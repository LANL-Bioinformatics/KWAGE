#include <mpi.h>

#ifdef _OPENMP
#include <omp.h> // For omp_get_max_threads()
#endif // _OPENMP

#include <math.h>

#include <sstream>
#include <fstream>
#include <unordered_set>

// DEBUG
#include <iostream>

#include "database.h"
#include "file_util.h"
#include "bigsi++.h"
#include "sort.h"

#include "slice_z.h"

using namespace std;

// Global variables for MPI
extern int mpi_numtasks;
extern int mpi_rank;

// From options.cpp (the allowed list of Bloom filter database file extensions)
extern const char* allowed_db_extentions[];

// For MAX_LOC MPI reduction
struct FloatInt{
	float value;
	int rank;
};

void merge_bloom_filters(const deque<BloomFilter> &m_db, const BloomParam &m_param, 
	UpdateInfo &m_progress,
	const BuildOptions &m_opt)
{
	time_t profile = time(NULL);
	
	stringstream ssout;
	
	ssout << m_opt.output_dir << PATH_SEPARATOR << m_param.kmer_len << "_" << hash_name(m_param.hash_func);
	
	if(mpi_rank == 0){
	
		if( !make_dir( ssout.str() ) ){
			
			cerr << "Error creating: " << ssout.str() << endl;
			throw __FILE__ ":merge_bloom_filters: Unable to create output directory (0)";
		}
	}
	
	ssout << PATH_SEPARATOR << m_param.log_2_filter_len;
	
	if(mpi_rank == 0){
	
		if( !make_dir( ssout.str() ) ){
			
			cerr << "Error creating: " << ssout.str() << endl;
			throw __FILE__ ":merge_bloom_filters: Unable to create output directory (1)";
		}
	}
	
	ssout << PATH_SEPARATOR << m_param.num_hash;
	
	if(mpi_rank == 0){
	
		if( !make_dir( ssout.str() ) ){
			
			cerr << "Error creating: " << ssout.str() << endl;
			throw __FILE__ ":merge_bloom_filters: Unable to create output directory (2)";
		}
	}
	
	// Count the number of filters (amoung all ranks)
	unsigned int count = m_db.size();
	
	// This MPI_Allreduce function acts as a barrier to allow the rank 0 process to 
	// construct any required output directories before the rank > 0 need to navigate
	// these directories.
	MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
	
	if(count == 0){
		
		// There are no bloom filters to write, something has gone wrong!
		// (there should always be at least one filter to write).
		throw __FILE__ ":merge_bloom_filters: Did not find any Bloom filters to write!";
	}
	
	const string target_dir = ssout.str();
	
	deque< pair<string, unsigned int> > existing_files;
	
	if(mpi_rank == 0){
		
		// Rank 0 reads any existing data files and counts the number of Bloom filters in each file
		// (to avoid having multiple ranks access the same files at the same time).
		FindFiles ff(target_dir);

		while( ff.next() ){

			if( find_file_extension(ff.name(), allowed_db_extentions) ) {
				
				// Read the number of Bloom filters for this file
				ifstream fin(ff.name().c_str(), ios::binary);
					
				if(!fin){
					throw __FILE__ ":merge_bloom_filters: Unable to open Bloom filter file";
				}
					
				DBFileHeader header;
				
				binary_read(fin, header);
					
				if(!fin){
					throw __FILE__":merge_bloom_filters: Error reading database file header";
				}
					
				existing_files.push_back( make_pair(ff.name(), header.num_filter) );
			}
		}
	}
	
	broadcast(existing_files, mpi_rank, 0);
	
	const size_t num_existing_files = existing_files.size();
	
	m_progress << "\tFound " << num_existing_files << " existing file" 
		<< ( (num_existing_files != 1) ? "s" : "" );
	m_progress.flush();
	m_progress.close();
	
	// The new total number of Bloom filters that must be stored
	unsigned int num_filter = count;
	
	for(deque< pair<string, unsigned int> >::const_iterator i = existing_files.begin();
		i != existing_files.end();++i){
		
		num_filter += i->second;
	}
	
	m_progress << "\tFound " << num_filter << " Bloom filter" 
		<< ( (num_filter != 1) ? "s" : "" );
	m_progress.flush();
	m_progress.close();
	
	bool use_compression = m_opt.compress && (num_filter >= MIN_NUM_FILTERS_FOR_COMPRESSION);
	
	// Represent each Bloom filter with a "partial" filter that only stores the first N bits.
	// Testing using Bacillus cereus group genomes suggests that the first (approximately) 320 bits
	// of a filter is sufficient to compute accurate distances between filters (the linear correlation 
	// coefficient between the partial filter distances and the full filter distances is approximately 0.96).
	
	const size_t filter_len = m_param.filter_len();
	
	deque<SubFilter> partial;

	// Create subfilters from the Bloom filters stored in the memory of this MPI rank
	for(deque<BloomFilter>::const_iterator i = m_db.begin();i != m_db.end();++i){

		partial.push_back( SubFilter() );

		partial.back().assign_bits(m_opt.subfilter_bits, *i);
		partial.back().set_src_loc( i - m_db.begin() );
	}

	m_progress << "\tLoading partial filters";
	m_progress.flush();
	
	#ifdef _OPENMP
		const int max_num_threads = omp_get_max_threads();
	#else
		const int max_num_threads = 1;
	#endif // _OPENMP
		
	//////////////////////////////////////////////////////////////////////////////////////
	// Initialize the compression and decompression engines (regardless if we are
	// using compression -- for now). Since the compression (and decompression) engines
	// are not thread safe (and there is some initialization overhead), we will create
	// a compression engine for every possible thread now.
	
	// Since the decompression step is currently closely tied to reading data from disk,
	// this step is still serial. If and when we switch to memory mapping database files for
	// fast reading, we can perform slice decompression in parallel.
	
	InflateSlice<MAX_COMPRESSED_BYTES> decompressor;
	
	vector< CompressSlice<MAX_COMPRESSED_BYTES> > compressor(max_num_threads);
	
	// Track the efficiency of our compression scheme
	size_t final_uncompressed_bytes = 0;
	size_t final_compressed_bytes = 0;

	// In addition to the in memory Bloom filters, each rank is also responsible for 
	// creating subfilters for a unique subset of database files
	for(size_t i = 0;i < num_existing_files;++i){
		
		// Skip the files assigned to other ranks
		if( (int)(i%mpi_numtasks) != mpi_rank){
			continue;
		}
		
		ifstream fin(existing_files[i].first.c_str(), ios::binary);

		if(!fin){
			throw __FILE__ ":merge_bloom_filters: Unable to open file";
		}

		DBFileHeader header;

		binary_read(fin, header);

		if(!fin){
			throw __FILE__":merge_bloom_filters: Error reading database file header (2)";
		}
		
		if(header.log_2_filter_len != m_param.log_2_filter_len){
			throw __FILE__ ":merge_bloom_filters: Existing file has an unexpected filter length";
		}
		
		if(header.num_hash != m_param.num_hash){
			throw __FILE__ ":merge_bloom_filters: Existing file has an unexpected number of hash functions";
		}
				
		if(filter_len < m_opt.subfilter_bits){
			throw __FILE__ ":merge_bloom_filters: bloom_filter_len < subfilter_bits";
		}
		
		// Allocate space to store the partial filters in this file
		const size_t partial_start = partial.size();

		for(uint32_t j = 0;j < header.num_filter;++j){

			partial.push_back( SubFilter() );

			SubFilter &ref = partial.back();

			ref.resize(m_opt.subfilter_bits);
			ref.unset_all_bits(); // Clear all bits
			ref.set_file(i);
			ref.set_src_loc(j);
		}

		BitVector slice(header.num_filter);

		if(header.compression != NO_COMPRESSION){

			// Read the compressed size of all bitslices in the file
			unsigned char *bitslice_compressed_len = new unsigned char[filter_len];

			fin.read( (char*)bitslice_compressed_len, filter_len );

			if(!fin){
				throw __FILE__":merge_bloom_filters: Error reading compressed slice length";
			}
				
			for(size_t j = 0;j < m_opt.subfilter_bits;++j){
				
				if(bitslice_compressed_len[j] == 0){
					
					// Compression failed to generate a smaller bit slice, so
					// we reverted to writing an uncompressed bit slice
					slice.read(fin);
					
					if(!fin){
						throw __FILE__":merge_bloom_filters: Error reading slice from file";
					}
				}
				else{

					decompressor.inflate(fin, bitslice_compressed_len[j]);

					slice.read( decompressor.ptr() );
				}

				// Pack this bit slice
				for(uint32_t k = 0;k < header.num_filter;++k){

					// Since the bits in partial list have all been
					// initialized to false, we only need to set the
					// bits that are true.
					if( slice.get_bit(k) ){
						partial[partial_start + k].set_bit(j);
					}
				}
			}

			if(bitslice_compressed_len != NULL){

				delete [] bitslice_compressed_len;
				bitslice_compressed_len = NULL;
			}
		}
		else{

			// If the bitslices are not compressed, then the size of each bitslice
			// is the length (in bytes) of a BitVector with header.num_filter bits
			for(size_t j = 0;j < m_opt.subfilter_bits;++j){				
				
				slice.read(fin);

				if(!fin){
					throw __FILE__":merge_bloom_filters: Error reading slice from file (2)";
				}
					
				// Pack this bit slice
				for(uint32_t k = 0;k < header.num_filter;++k){

					// Since the bits in partial list have all been
					// initialized to false, we only need to set the
					// bits that are true.
					if( slice.get_bit(k) ){
						partial[partial_start + k].set_bit(j);
					}
				}
			}
		}
	}

	m_progress << "\tComputing Bloom filter centroid";
	m_progress.flush();
	
	// Compute the center of the distribution of the partial Bloom filters.
	// The center is the average value of each bit
	float *ave = new float [m_opt.subfilter_bits];

	if(ave == NULL){
		throw __FILE__ ":merge_bloom_filters: Unable to allocate ave";
	}

	memset( ave, 0, m_opt.subfilter_bits*sizeof(float) );
	
	for(deque<SubFilter>::const_iterator i = partial.begin();i != partial.end();++i){
	
		for(uint32_t j = 0;j < m_opt.subfilter_bits;++j){

			if( i->get_bit(j) ){
				++ave[j];
			}
		}
	}

	// Share bit count values with all ranks
	unsigned int bit_norm = partial.size();

	MPI_Allreduce(MPI_IN_PLACE, &bit_norm, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, ave, m_opt.subfilter_bits, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	if(bit_norm > 0){
		
		const float norm = 1.0f/bit_norm;

		for(uint32_t i = 0;i < m_opt.subfilter_bits;++i){
			ave[i] *= norm;
		}
	}

	// Find the single partial filter that is farthest from the average
	SubFilter edge_filter;
	float edge_distance2 = 0.0;
	const size_t num_partial = partial.size();

	for(size_t i = 0;i < num_partial;++i){

		float d2 = 0.0f;

		for(uint32_t j = 0;j < m_opt.subfilter_bits;++j){
			
			const float delta = partial[i].get_bit(j) - ave[j];

			d2 += delta*delta;
		}

		// Use ">=" for the special case of a file with a single Bloom filter
		// or multiple identical Bloom filters
		if(d2 >= edge_distance2){
			
			edge_distance2 = d2;
			edge_filter = partial[i];
		}
	}

	// Since we have identified the filter that is farthest from the center of
	// all filters, we can delete the average.
	if(ave != NULL){

		delete [] ave;
		ave = NULL;
	}

	// Find the most extreme partial filter across all ranks
	FloatInt max_loc;

	max_loc.value = edge_distance2;
	max_loc.rank = mpi_rank;

	MPI_Allreduce(MPI_IN_PLACE, &max_loc, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD);

	m_progress << "\tFound most extreme Bloom filter " 
		<< sqrt(max_loc.value) << " from the centroid";
	m_progress.flush();
	
	// The rank with the most extreme partial filter needs to share it with
	// the remaining ranks
	broadcast(edge_filter, mpi_rank, max_loc.rank);

	// To order all filters, we find the filter that minimizes the distance
	// to both the edge_filter and the last filter added:
	// weighted distance = gamma*(distance to edge) + (1 - gamma)*(distance to last)
	//const float gamma = 0.5f; // The effect of this weight on search and compression must be tested!!!!

	// |Gamma|Bacillus DB size|
	// 1.0  4.0G
	// 0.8  
	// 0.7  
	// 0.6  2.9G 
	// 0.5	2.9G
	// 0.4  2.9G	
	// 0.3  2.8G
	// 0.2  2.8G
 	// 0.1	2.8G
 	// 0.0	2.8G
	// Rank every partial filter on all ranks based on the above mentioned criteria.
	uint32_t partial_order_index = 0;

	const float largest_distance = float(m_opt.subfilter_bits);

	m_progress << "\tOrdering Bloom filters by similarity";
	m_progress.flush();
	
	while(true){

		// The max distance between any two filters is
		// m_opt.subfilter_bits, so initialize the distance to
		// a greater value
		float best_distance = largest_distance + 1.0f;

		size_t best_partial_index = 0;

		for(size_t i = 0;i < num_partial;++i){

			// Does the current partial filter still need to be ranked?
			if(partial[i].get_dst_loc() == SubFilter::INVALID_LOC){
			
				float d2 = 0.0f;

				for(uint32_t j = 0;j < m_opt.subfilter_bits;++j){
				
					float delta = partial[i].get_bit(j) - edge_filter.get_bit(j);
					
					d2 += delta*delta;
				}

				if(best_distance > d2){
					
					best_distance = d2;
					best_partial_index = i;
				}
			}
		}

		FloatInt min_loc;

		min_loc.value = best_distance;
		min_loc.rank = mpi_rank;

		MPI_Allreduce(MPI_IN_PLACE, &min_loc, 1, MPI_FLOAT_INT, MPI_MINLOC, MPI_COMM_WORLD);

		if(min_loc.value > largest_distance){
			break;
		}
		
		if(min_loc.rank == mpi_rank){
			
			partial[best_partial_index].set_dst_loc(partial_order_index);
			edge_filter = partial[best_partial_index];

			// DEBUG
			//cout << partial_order_index << '\t' << min_loc.value << '\t'
			//	<< m_db[ partial[best_partial_index].get_src_loc() ].get_info() << endl;
		}

		// The new edge filter selected must be shared with all ranks	
		broadcast(edge_filter, mpi_rank, min_loc.rank);

		++partial_order_index;
	}
	
	// It is now, finally, time to write the Bloom filter files to disk
	// We won't delete the old files (if any) until we have confirmed the
	// integrity of the new files. Rank zero will be responsible for collecting, optionally
	// compressing and then writing the Boom filters to disk.

	// Order the partial filters by their destination location
	SORT( partial.begin(), partial.end() );
	
	// Before writing the new database files, free up as much space as possible
	// by deallocating the bitslices associated with all of the partial
	// filters (but keep the source and destination information for each
	// partial filter)
	for(deque<SubFilter>::iterator i = partial.begin();i != partial.end();++i){
		i->clear(); // This only frees the bitvector buffer
	}

	// Enable fast mapping between the global filter order and the local order on this rank
	unordered_map<unsigned int, unsigned int> global_to_local_filter;
	
	for(unsigned int i = 0;i < num_partial;++i){
		global_to_local_filter[partial[i].get_dst_loc()] = i;
	}
	
	unsigned int curr_file_index = 0;

	// Since the maximum number of Bloom filters per file is NUM_FILTER_CHUNK, we will be
	// creating a new file for each iteration of this for-loop.
	for(unsigned int i = 0;i < num_filter;i += NUM_FILTER_CHUNK){

		++curr_file_index;
	
		// Each output database file contains at most NUM_FILTER_CHUNK Bloom filters
		const unsigned int filters_per_file = min(num_filter - i, NUM_FILTER_CHUNK);

		// The first filter index
		const unsigned int filter_begin = i;

		// One past the last filter index
		const unsigned int filter_end = i + filters_per_file;

		// Allocate a buffer for collecting the bit-slice data for all Bloom filters from
		// each of the ranks. To do this, we need to know the number of bytes in each
		// bit slice
		const unsigned int bytes_per_bitslice = filters_per_file/BitVector::BITS_PER_BLOCK +
			( (filters_per_file%BitVector::BITS_PER_BLOCK) > 0 ? 1 : 0);

		// Since we always choose Bloom filter lengths that are a power of two,
		// make the number of slices per buffer a power of two.
		const unsigned int num_buffer_slice = 2048;

		unsigned char* buffer = new unsigned char[num_buffer_slice * bytes_per_bitslice];

		if(buffer == NULL){
			throw __FILE__ ":merge_bloom_filters: Unable to allocate bit-slice collection buffer";
		}

		unsigned char *compressed_bitslice_len = NULL;
		
		if(use_compression){
				
			compressed_bitslice_len = new unsigned char [filter_len];
				
			if(compressed_bitslice_len == NULL){
				throw __FILE__ ":merge_bloom_filters: Unable to allocate bitslice length place holder array";
			}
				
			memset( compressed_bitslice_len, 0, filter_len*sizeof(unsigned char) );
		}

		// For now, only rank 0 will actually write the database files
		ofstream fout;

		stringstream ssout;

		// The output filename
		ssout << target_dir << PATH_SEPARATOR << "new_" << curr_file_index << ".db";
		
		DBFileHeader output_header;
		
		output_header.kmer_len = m_param.kmer_len;
		output_header.num_hash = m_param.num_hash;
		output_header.log_2_filter_len = m_param.log_2_filter_len;
		output_header.num_filter = filters_per_file;
		output_header.compression = use_compression ? ZLIB_RLE_COMPRESSION : NO_COMPRESSION;

		if(mpi_rank == 0){

			fout.open(ssout.str().c_str(), ios::binary);

			if(!fout){
				throw __FILE__ ":merge_bloom_filters: Unable to open bloom filter for writing";
			}

			// Write the header information for this file
			binary_write(fout, output_header);

			if(!fout){
				throw __FILE__ ":merge_bloom_filters: Error writing database file header";
			}
			
			// If we are compressing the bitslices, we will need to make space in the output
			// file to store the compressed lengths of each bitslice.
			if(use_compression){
				
				fout.write( (char*)compressed_bitslice_len, filter_len*sizeof(unsigned char) );
				
				if(!fout){
					throw __FILE__ ":merge_bloom_filters: Error writing compressed slice length placeholder";
				}
			}
		}

		// Prepare to read Bloom filter data from all of the existing files that this rank
		// is responsible for. Note that older C++ compilers do not allow us to create a 
		// vector of ifstream objects (due to the lack of a copy constructor for these objects).
		// To increase our portability, we will fall back to C-style FILE* pointers.
		
		vector<ifstream*> fin_ptr(num_existing_files, NULL);
		vector<uint32_t> filters_per_target_file(num_existing_files);
		vector<unsigned char*> existing_compressed_slice_len(num_existing_files, NULL);
		
		unordered_set<unsigned int> target_files;
		
		for(unsigned int j = filter_begin;j < filter_end;++j){			
			
			unordered_map<unsigned int, unsigned int>::const_iterator iter = global_to_local_filter.find(j);
			
			if( (iter != global_to_local_filter.end() ) && (partial[iter->second].get_file() != SubFilter::INVALID_FILE) ){
			    	target_files.insert( partial[iter->second].get_file() );
			}
		}

		const unsigned int num_slice_chunk = filter_len/num_buffer_slice;
	
		if( (num_slice_chunk == 0) || (filter_len%num_buffer_slice != 0) ){
			throw __FILE__ ":merge_bloom_filters: Unexpected remainder (non-power of two) or too small number of bit slices in file";
		}
		
		const unsigned int update_slice_every = max(1U, num_slice_chunk/100);
		bool disable_compression = false;
		
		for(unsigned int j = 0;j < num_slice_chunk;++j){
	
			if(disable_compression){
				
				use_compression = false;
				disable_compression = false;
				
				// Start all over again *without* compression
				j = 0;
				
				// Modify the output header to show that we have turned off compression
				output_header.compression = use_compression ? ZLIB_RLE_COMPRESSION : NO_COMPRESSION;
				
				if(compressed_bitslice_len != NULL){
				
					delete [] compressed_bitslice_len;
					compressed_bitslice_len = NULL;
				}

				if(mpi_rank == 0){
					
					// Reopen the output file and rewrite the header
					fout.close();
					
					fout.open(ssout.str().c_str(), ios::binary);

					if(!fout){
						throw __FILE__ ":merge_bloom_filters: Unable to open bloom filter for writing (2)";
					}

					// Write the header information for this file
					binary_write(fout, output_header);

					if(!fout){
						throw __FILE__ ":merge_bloom_filters: Error writing database file header (2)";
					}
				}
			}
			
			// Open the input files on the first iteration, or if we disable compression and restart
			if(j == 0){
			
				for(unordered_set<unsigned int>::const_iterator k = target_files.begin();
					k != target_files.end();++k){

					if(fin_ptr[*k] != NULL){
					
						fin_ptr[*k]->close();
						fin_ptr[*k] = NULL;
					}

					fin_ptr[*k] = new ifstream;

					if(fin_ptr[*k] == NULL){
						throw __FILE__ ":merge_bloom_filters: Unable to allow ifstream";
					}

					ifstream &fin_ref = *(fin_ptr[*k]);

					fin_ref.open(existing_files[*k].first.c_str(), ios::binary);

					if(!fin_ref){
						throw __FILE__ ":merge_bloom_filters: Unable to open existing file for reading Bloom filters";
					}

					// Read the header
					DBFileHeader header;

					binary_read(fin_ref, header);

					if(!fin_ref){
						throw __FILE__":merge_bloom_filters: Error reading database file header (3)";
					}

					if(header.log_2_filter_len != m_param.log_2_filter_len){
						throw __FILE__ ":merge_bloom_filters: Unexpected mismatch in Bloom filter lengths";
					}

					filters_per_target_file[*k] = header.num_filter;

					if(existing_compressed_slice_len[*k] != NULL){
					
						delete [] existing_compressed_slice_len[*k];
						existing_compressed_slice_len[*k] = NULL;
					}

					if(header.compression != NO_COMPRESSION){

						// As this is a compressed file, we will need to store the compressed lengths
						// of each bitslice
						existing_compressed_slice_len[*k] = new unsigned char [filter_len];

						if(existing_compressed_slice_len[*k] == NULL){
							throw __FILE__ ":merge_bloom_filters: Unable to allocate compressed slice length for existing file";
						}

						fin_ref.read( (char*)existing_compressed_slice_len[*k], filter_len );

						if(!fin_ref){
							throw __FILE__":merge_bloom_filters: Error reading compressed slice lengths (2)";
						}
					}
				}
			}
			
			if(j%update_slice_every == 0){
			
				m_progress << "\tWriting " << (use_compression ? "compressed " : "") 
					<< "bitslice data for file " << curr_file_index << " " 
					<< (100.0*j)/num_slice_chunk << "%";
				m_progress.flush();
			}
			
			const unsigned int slice_begin = j*num_buffer_slice;
			const unsigned int slice_end = slice_begin + num_buffer_slice;
			
			// Initialize the buffer to all zeros
			memset(buffer, 0, num_buffer_slice * bytes_per_bitslice);
			
			// Write the Bloom filter slices that are currently in memory to the buffer
			// (each rank performs this task independently).
			for(unsigned int g = filter_begin;g < filter_end;++g){
				
				unordered_map<unsigned int, unsigned int>::const_iterator iter = global_to_local_filter.find(g);
				
				if( iter == global_to_local_filter.end() ){
					continue;
				}
				
				if(partial[iter->second].get_file() == SubFilter::INVALID_FILE){
					
					// This refers to an in-memory Bloom filter
					const BloomFilter &ref = m_db[ partial[iter->second].get_src_loc() ];
					
					const unsigned int delta_g = g - filter_begin;
					
					// Write (in transposed orientation) the filter elements
					// [slice_begin, slice_end) to the buffer for this group
					unsigned char* ptr = buffer + delta_g/BitVector::BITS_PER_BLOCK;
					const unsigned int offset = delta_g%BitVector::BITS_PER_BLOCK;
					
					for(unsigned int k = slice_begin;k < slice_end;++k, ptr += bytes_per_bitslice){
						*ptr |= ref.get_bit(k) << offset;
					}
				}
			}
			
			// Now, for any Bloom filter whose source is an existing file, we need to load
			// and then pack the bits (as with the slices in memory, each rank can perform
			// this task independently).
			for(unordered_set<unsigned int>::const_iterator f = target_files.begin();f != target_files.end();++f){
				
				ifstream &fin_ref = *(fin_ptr[*f]);
				
				BitVector slice(filters_per_target_file[*f]);
				
				unsigned char* ptr = buffer;
				
				for(unsigned int k = slice_begin;k < slice_end;++k, ptr += bytes_per_bitslice){
				
					if( (existing_compressed_slice_len[*f] == NULL) ||
					    (existing_compressed_slice_len[*f][k] == 0) ){ // Uncompressed

						slice.read(fin_ref);

						if(!fin_ref){
							throw __FILE__ ":merge_bloom_filters: Unable to read uncompressed slice";
						}
					}
					else{ // Compressed
						
						// Read the compressed data from the input file
						decompressor.inflate(fin_ref, existing_compressed_slice_len[*f][k]);

						slice.read( decompressor.ptr() );
					}

					for(unsigned int g = filter_begin;g < filter_end;++g){

						unordered_map<unsigned int, unsigned int>::const_iterator iter = global_to_local_filter.find(g);

						if( iter == global_to_local_filter.end() ){
							continue;
						}
						
						const SubFilter &ref = partial[iter->second];
						
						if(ref.get_file() == *f){

							// This group is found in the current file. Since we are reading the
							// slice directly from a database file, we do *not* need to transpose it!
							
							const unsigned int delta_g = g - filter_begin;
							
							unsigned char* ptr2 = ptr + delta_g/BitVector::BITS_PER_BLOCK;
							const unsigned int offset = delta_g%BitVector::BITS_PER_BLOCK;

							*ptr2 |= slice.get_bit( ref.get_src_loc() ) << offset;
						}
					}
				}
			}

			if(use_compression){

				// Share the current slice buffer amoung all ranks
				MPI_Allreduce(MPI_IN_PLACE, buffer, num_buffer_slice * bytes_per_bitslice, 
					MPI_BYTE, MPI_BOR, MPI_COMM_WORLD);

				// Each rank is responsible for compressing a subset of slices. Each rank is
				// assigned a unique set of slices, which can be compressed in parallel
				// by multiple threads on each rank.
				//
				// Note that for this compression process, it is essential to keep the memory
				// access as local as possible! Assign each MPI rank a *continguous* block of 
				// slices to compress. (The initial version that used non-local memory access
				// was super-super slow!)
				
				const unsigned int chunk = max(1U, num_buffer_slice/mpi_numtasks);
				
				const unsigned int local_begin = chunk*mpi_rank;
				const unsigned int local_end = ( mpi_rank == (mpi_numtasks - 1) ) ?
					num_buffer_slice : // Remainder slices are assigned to the last rank
					local_begin + chunk;
				
				// Zero the bits in the buffer that this rank is *not* going to compress to
				// allow the final buffer to be combined with a binary OR. 
				memset(buffer, 0, bytes_per_bitslice*local_begin); // Zero the bits *before* local_begin
				
				// Zero the buffer bits *on or after* local_end.
				memset( buffer + bytes_per_bitslice*local_end, 0, 
					bytes_per_bitslice*(num_buffer_slice - local_end) ); 
				
				#pragma omp parallel
				{
					#ifdef _OPENMP
						const int tid = omp_get_thread_num();
					#else
						const int tid = 0;
					#endif // _OPENMP

					CompressSlice<MAX_COMPRESSED_BYTES> &compressor_ref = compressor[tid];
					
					#pragma omp for					
					for(unsigned int k = local_begin;k < local_end;++k){
						
						// All buffer access is relative to the slice_begin index
						unsigned char* ptr = buffer + bytes_per_bitslice*k;

						// This rank is responsible for compressing the current slice
						if( compressor_ref.compress(ptr, bytes_per_bitslice) ){

							// Access into the compressed_bitslice_len array is
							// relative to slice_begin.
							const unsigned int index = k + slice_begin;
							
							// We successfully compressed the buffer smaller than the input size of
							// bytes_per_bitslice. Write this compressed slice to disk
							compressed_bitslice_len[index] = compressor_ref.size();

							memcpy(ptr, compressor_ref.ptr(), compressed_bitslice_len[index]);
						}
					}
				}
				
				// Copy the compressed slice lengths to rank 0
				if(mpi_rank == 0){	
					MPI_Reduce(MPI_IN_PLACE, compressed_bitslice_len + slice_begin, slice_end - slice_begin, 
						MPI_BYTE, MPI_MAX, 0, MPI_COMM_WORLD);
				}
				else{
					MPI_Reduce(compressed_bitslice_len + slice_begin, NULL, slice_end - slice_begin, 
						MPI_BYTE, MPI_MAX, 0, MPI_COMM_WORLD);
				}
			}

			if(mpi_rank == 0){
				
				// Merge the buffer data from all ranks using the binary OR operator
				MPI_Reduce(MPI_IN_PLACE, buffer, num_buffer_slice * bytes_per_bitslice, 
					MPI_BYTE, MPI_BOR, 0, MPI_COMM_WORLD);
	
				if(use_compression){
					
					// Attempt to compress each row (i.e. slice) of the buffer. If the compression does
					// not result in a smaller number of bytes, then write the uncompressed buffer row
					// instead
					unsigned char* ptr = buffer;
				
					for(unsigned int k = slice_begin;k < slice_end;++k, ptr += bytes_per_bitslice){
						
						final_uncompressed_bytes += bytes_per_bitslice;

						// Did we successfully compress the bit slice?
						if(compressed_bitslice_len[k] == 0){
							
							// We were unable to efficiently compress the input bitslice, fall back to 
							// writing the uncompressed data
							fout.write( (char*)ptr, bytes_per_bitslice*sizeof(unsigned char) );

							final_compressed_bytes += bytes_per_bitslice;
						}
						else{

							// We successfully compressed the buffer smaller than the input size of
							// bytes_per_bitslice. Write this compressed slice to disk
							fout.write( (char*)ptr, compressed_bitslice_len[k]*sizeof(unsigned char) );

							final_compressed_bytes += compressed_bitslice_len[k];
						}

						if(!fout){
							throw __FILE__ ":merge_bloom_filters: Error writing slice to disk";
						}
					}
				}
				else{

					fout.write( (char*)buffer, num_buffer_slice*bytes_per_bitslice*sizeof(unsigned char) );
					
					if(!fout){
						throw __FILE__ ":merge_bloom_filters: Error writing slice buffer to disk";
					}
				}
			}
			else{

				// Merge the buffer data from all ranks to rank 0 using the binary OR operator
				MPI_Reduce(buffer, NULL, num_buffer_slice * bytes_per_bitslice, MPI_BYTE, MPI_BOR, 0, MPI_COMM_WORLD);
			}
			
			if(use_compression){
				
				float compression_level = 0.0;
				
				// The amount of compression that we can obtain depends on the kmer-composition 
				// of the input Bloom filters. To avoid wasting time compressing Bloom filter slices that
				// do not compress well, perform a quick check the compression level obtained for the 
				// partial Bloom filters. If we don't obtain the desired compression level (or better),
				// disable compression (which will dramatically speed up the database contruction step!).
				if(mpi_rank == 0){
                                        compression_level = double(final_compressed_bytes)/final_uncompressed_bytes;
                                }
			
				// Only rank 0 tracks the final_compressed_bytes and final_uncompressed_bytes
				// needed to compute the compression level	
				MPI_Bcast(&compression_level, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

				if(compression_level > m_opt.compression_threshold){
					
					disable_compression = true;
					
					m_progress << "\tDisabling compression (compression level " 
						<< compression_level << " > threshold of " 
						<< m_opt.compression_threshold << ")";
					m_progress.flush();
					m_progress.close();
				}
			}
		}
		
		// We are finished writing the bitslices for this file (yay!).
		if(use_compression){

			m_progress << "\tCompression reduced size to " 
				<< (100.0*final_compressed_bytes)/final_uncompressed_bytes << "% of original";
			m_progress.flush();
			m_progress.close();
		}

		delete [] buffer;
		buffer = NULL;
				
		for(vector<unsigned char*>::iterator j = existing_compressed_slice_len.begin();
			j != existing_compressed_slice_len.end();++j){
			
			if(*j != NULL){
				
				delete [] *j;
				*j = NULL;
			}
		}
		
		m_progress << "\tWriting metadata for file " << curr_file_index;
		m_progress.flush();
			
		// We now need to write the group information data (i.e. the positions of each group info string in the
		// file, and the corresponding string data). To do this, rank 0 (which is responsible for the file writing)
		// writes a *dummy* array to make room for the group location information.
		unsigned long int *group_loc_buffer = NULL;
		
		// Record the output file location of the group loc data
		unsigned long int fout_loc_group_loc = 0;
		
		// Make a hole that will be later filled with the
		// file locations of the filter information strings
		if(mpi_rank == 0){
			
			fout_loc_group_loc = fout.tellp();

			group_loc_buffer = new unsigned long int [filters_per_file];
			
			if(group_loc_buffer == NULL){
				throw __FILE__ ":merge_bloom_filters: Unable to allocate group_loc_buffer";
			}
			
			memset( group_loc_buffer, 0, filters_per_file*sizeof(unsigned long int) );
			
			fout.write( (char*)group_loc_buffer, filters_per_file*sizeof(unsigned long int) );
			
			if(!fout){
				throw __FILE__ ":merge_bloom_filters: Error writing dummy group loc buffer";
			}
		}
		
		// Read the input file locations for each of the input file-specific group info locations
		vector<unsigned long int*> input_group_loc(num_existing_files, NULL);
		
		for(unordered_set<unsigned int>::const_iterator j = target_files.begin();
			j != target_files.end();++j){
			
			input_group_loc[*j] = new unsigned long int[ filters_per_target_file[*j] ];
			
			if(input_group_loc[*j] == NULL){
				throw __FILE__ ":merge_bloom_filters: Unable to allocate group information string location array";
			}
			
			ifstream &fin_ref = *(fin_ptr[*j]);
			
			fin_ref.read( (char*)input_group_loc[*j], filters_per_target_file[*j]*sizeof(unsigned long int) );
			
			if(!fin_ref){
				throw __FILE__ ":merge_bloom_filters: Unable to read group information string location";
			}
		}
		
		for(unsigned int j = filter_begin;j < filter_end;++j){
		
			// To prevent the race condition of ranks returning "out of order" filter info
			// strings, force all ranks to proceed at the pace of rank 0
			broadcast(j, mpi_rank, 0);

			// Which ever rank "owns" this group information needs to send it to rank 0 for
			// writing to disk
			unordered_map<unsigned int, unsigned int>::const_iterator iter = global_to_local_filter.find(j);
			
			FilterInfo info;
			
			if(mpi_rank == 0){
				
				if( iter == global_to_local_filter.end() ){
					
					// Another rank has this info metadata
					unsigned int info_buffer_len;
			
					MPI_Status status;
					
					if(MPI_Recv(&info_buffer_len, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, INFO_LEN, MPI_COMM_WORLD, &status) != MPI_SUCCESS){
						throw __FILE__ ":merge_bloom_filters: Error receiving INFO_STRING_LEN msg";
					}

					unsigned char *buffer = new unsigned char[info_buffer_len];

					if(buffer == NULL){
						throw __FILE__ ":merge_bloom_filters: Unable to allocate info string buffer";
					}

					if(MPI_Recv(buffer, info_buffer_len, MPI_BYTE, status.MPI_SOURCE, INFO_BUFFER, MPI_COMM_WORLD, &status) != MPI_SUCCESS){
						throw __FILE__ ":merge_bloom_filters: Error receiving info string data";
					}
					
					mpi_unpack(buffer, info);
					
					delete [] buffer;
					buffer = NULL;
				}
				else{

					const unsigned int file_index = partial[iter->second].get_file();
					
					if(file_index == SubFilter::INVALID_FILE){
					
						// In memory (on rank 0) Bloom filter
						info = m_db[partial[iter->second].get_src_loc()].get_info();
					}
					else{
						// Load from disk
						ifstream &fin_ref = *(fin_ptr[file_index]);
						
						fin_ref.seekg(input_group_loc[file_index][partial[iter->second].get_src_loc()]);

						binary_read(fin_ref, info);
					}
				}
				
				// Group numbering starts at zero within each file
				group_loc_buffer[j - filter_begin] = fout.tellp();
				
				binary_write(fout, info);

				if(!fout){
					throw __FILE__ ":merge_bloom_filters: Error writing filter information";
				}
			}
			else{

				if( iter == global_to_local_filter.end() ){
					
					// A different rank has this filter metadata
					continue;
				}
				
				const unsigned int file_index = partial[iter->second].get_file();
				
				if(file_index == SubFilter::INVALID_FILE){
					
					// In memory Bloom filter
					info = m_db[partial[iter->second].get_src_loc()].get_info();
				}
				else{
					// Load from disk
					ifstream &fin_ref = *(fin_ptr[file_index]);
					
					fin_ref.seekg(input_group_loc[file_index][partial[iter->second].get_src_loc()]);
					
					binary_read(fin_ref, info);
				}
				
				// Send this filter metadata to rank 0
				const unsigned int buffer_len = mpi_size(info);
				
				if(MPI_Send( (void*)&buffer_len, 1, MPI_UNSIGNED, 0, INFO_LEN, MPI_COMM_WORLD ) != MPI_SUCCESS){
					throw __FILE__ ":merge_bloom_filters: Error sending msg length";
				}
				
				unsigned char* buffer = new unsigned char [buffer_len];
				
				if(buffer == NULL){
					throw __FILE__ ":merge_bloom_filters: Error allocating send buffer";
				}

				mpi_pack(buffer, info);
				
				if(MPI_Send( (void*)buffer, buffer_len, MPI_BYTE, 0, INFO_BUFFER, MPI_COMM_WORLD ) != MPI_SUCCESS){
					throw __FILE__ ":merge_bloom_filters: Error sending msg";
				}
		
				delete [] buffer;
				buffer = NULL;
			}
		}
		
		// Manually deallocate (which will close) the input files. This is required by our use of pointers to ifstream
		// objects (which is, in turn, required to retain compatibility with older C++ compilers which do not have a
		// copy constructor defined for ifstream).
		for(vector<ifstream*>::iterator j = fin_ptr.begin();j != fin_ptr.end();++j){
		
			if(*j != NULL){
			
				delete *j;
				*j = NULL;
			}
		}
		
		fin_ptr.clear();
		
		// The last step when writing a database file is to go back and update the locations of the group information strings
		if(mpi_rank == 0){
			
			fout.seekp(fout_loc_group_loc);
			
			fout.write( (char*)group_loc_buffer, filters_per_file*sizeof(unsigned long int) );
			
			if(!fout){
				throw __FILE__ ":merge_bloom_filters: Error writing actual group loc buffer";
			}
		}
		
		if(group_loc_buffer != NULL){
			
			delete [] group_loc_buffer;
			group_loc_buffer = NULL;
		}
		
		for(vector<unsigned long int*>::iterator j = input_group_loc.begin();j != input_group_loc.end();++j){
			
			if(*j != NULL){
			
				delete [] *j;
				*j = NULL;
			}
		}
		
		if(compressed_bitslice_len != NULL){
			
			// Write the compressed bitslice lengths (if needed)
			if(mpi_rank == 0){
				
				// Start writing immediately after the file header
				fout.seekp( sizeof(DBFileHeader) );
				
				fout.write( (char*)compressed_bitslice_len, filter_len*sizeof(unsigned char) );
				
				if(!fout){
					throw __FILE__ ":merge_bloom_filters: Error writing compressed slice length";
				}				
			}

			delete [] compressed_bitslice_len;
			compressed_bitslice_len = NULL;
		}
		
		m_progress << "\tDone writing file " << curr_file_index;
		m_progress.flush();
		m_progress.close();
	}
	
	if(mpi_rank == 0){
		
		if( !existing_files.empty() ){
		
			// Now we need to remove any existing files and rename the new files. 
			// Extra validation is needed!!!
			m_progress << "\tRemoving old database files";
			m_progress.flush();
			m_progress.close();
			
			for(deque< pair<string, unsigned int> >::const_iterator i = existing_files.begin();
				i != existing_files.end();++i){
				
				if(unlink( i->first.c_str() ) != 0){
					throw __FILE__ ":merge_bloom_filters: Error removing file!";
				}
			}
		}
		
		m_progress << "\tRenaming database files";
		m_progress.flush();
		m_progress.close();
		
		for(unsigned int i = 1;i <= curr_file_index;++i){
			
			stringstream ssold;
			stringstream ssnew;
			
			ssold << target_dir << PATH_SEPARATOR << "new_" << i << ".db";
			ssnew << target_dir << PATH_SEPARATOR << i << ".db";
			
			if(rename( ssold.str().c_str(), ssnew.str().c_str() ) != 0){
				throw __FILE__ ":merge_bloom_filters: Error renaming new file!";
			}
		}
		
		profile = time(NULL) - profile;
		
		m_progress << "\tDatabase complete in " << profile << " sec";
		m_progress.flush();
		m_progress.close();
	}
}
