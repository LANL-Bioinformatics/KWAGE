// Build a CALDERA Bloom filter database from input Bloom filter files
// J. D. Gans
// Bioscience Division, B-10
// Los Alamos National Laboratory
// Fri Oct 25 10:22:49 2019

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <zlib.h>

#include "maestro.h"
#include "caldera.h"
#include "mpi_util.h"
#include "file_util.h"
#include "binary_io.h"

using namespace std;

extern int mpi_rank;

bool build_db(const string &m_filename, const BloomParam &m_param, const deque<string> &m_bloom_files)
{
	try{
		
		const size_t num_filter = m_bloom_files.size();

		if(num_filter == 0){
			throw __FILE__ ":build_db: Empty Bloom filter inventory file";
		}

		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] Writing databaes file to " << m_filename << endl;
		cerr << "[" << mpi_rank << "] Found " << num_filter << " input Bloom filters" << endl;
		cerr << "[" << mpi_rank << "] log_2_filter_len = " << m_param.log_2_filter_len << endl;
		cerr << "[" << mpi_rank << "] num_hash = " << m_param.num_hash << endl;
		cerr << "[" << mpi_rank << "] hash_func = " << hash_name(m_param.hash_func) << endl;
		cerr << "[" << mpi_rank << "] kmer_len = " << m_param.kmer_len << endl;
		#endif // DEBUG_DB

		// Open all of the Bloom filter files and read the header information.
		// Note that older C++ compilers do not allow us to create a  vector of ifstream 
		// objects (due to the lack of a copy constructor for these objects).
		// To increase our portability, we will use *pointers* to ifstream objects
		// (which must be deleted!).
		vector<ifstream*> fin_ptr(num_filter, NULL);

		for(size_t i = 0;i < num_filter;++i){
			
			fin_ptr[i] = new ifstream();

			if(fin_ptr[i] == NULL){
				throw __FILE__ ":build_db: Unable to allocate ifstream for Bloom filter";
			}

			fin_ptr[i]->open(m_bloom_files[i].c_str(), ios::binary);

			if( !(*fin_ptr[i]) ){

				#ifdef DEBUG_DB
				cerr << "[" << mpi_rank << "] Unable to open Bloom filter: " << m_bloom_files[i] << endl;
				#endif // DEBUG_DB

				throw __FILE__ ":build_db: Unable to open Bloom filter file";
			}

			#ifdef DEBUG_DB
			cerr << "[" << mpi_rank << "] Reading Bloom filer from " << m_bloom_files[i] << endl;
			#endif // DEBUG_DB
		}

		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] All filter files have been opened." << endl;
		#endif // DEBUG_DB

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Make sure that all of the Bloom filters have been completely written
		for(size_t i = 0;i < num_filter;++i){

			unsigned char status;

			binary_read(*fin_ptr[i], status);

			if(status != BLOOM_MAGIC_COMPLETE){

				#ifdef DEBUG_DB
				cerr << "[" << mpi_rank << "] Incomplete Bloom filter: " << m_bloom_files[i] << endl;
				#endif // DEBUG_DB

				throw __FILE__ ":build_db: Incomplete Bloom filter";
			}
		}

		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] All Bloom filters are complete." << endl;
		#endif // DEBUG_DB

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read the BloomParam data from each filter. We won't be storing the this information for each filter,
		// but do need to make sure that the filter paramters are all the *same*.

		for(size_t i = 0;i < num_filter;++i){

			BloomParam local;

			binary_read(*fin_ptr[i], local);

			if( !(*fin_ptr[i]) ){

				#ifdef DEBUG_DB
				cerr << "[" << mpi_rank << "] Error reading BloomParam from " << m_bloom_files[i] << endl;
				#endif // DEBUG_DB

				throw __FILE__ ":build_db: Error reading BloomParam";
			}

			if(m_param != local){

				#ifdef DEBUG_DB
				cerr << "[" << mpi_rank << "] Inconsistent Bloom filter parameters for " << m_bloom_files[i] << endl;
				cerr << "[" << mpi_rank << "] kmer_len = " << local.kmer_len << "; expected " << m_param.kmer_len << endl;
				cerr << "[" << mpi_rank << "] log_2_filter_len = " << local.log_2_filter_len << "; expected " << m_param.log_2_filter_len << endl;
				cerr << "[" << mpi_rank << "] num_hash = " << local.num_hash << "; expected " << m_param.num_hash << endl;
				cerr << "[" << mpi_rank << "] hash_func = " << hash_name(local.hash_func) << "; expected " << hash_name(m_param.hash_func) << endl;
				#endif // DEBUG_DB

				throw __FILE__ ":build_db: Inconsistent Bloom parameters";
			}
		}

		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] All Bloom filter parameters are consistent." << endl;
		#endif // DEBUG_DB

		const size_t filter_len = m_param.filter_len();

		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read the CRC data from each filter. These crc32 values will be used to validate the Blooom filter
		// data that we read. After all Bloom filters have been read, the reference and computed CRC32 values
		// should be equal. Note that before we compute the running crc32 values, we must initialize
		// the crc32 value with a call to crc32_z(0L, Z_NULL, 0).
		vector< pair<unsigned int /*reference*/, unsigned int /*computed*/> > filter_crc32(num_filter, 
			make_pair(0, crc32_z(0L, Z_NULL, 0) ) );

		for(size_t i = 0;i < num_filter;++i){

			binary_read(*fin_ptr[i], filter_crc32[i].first);

			if( !(*fin_ptr[i]) ){

				#ifdef DEBUG_DB
				cerr << "[" << mpi_rank << "] Error reading CRC32 from " << m_bloom_files[i] << endl;
				#endif // DEBUG_DB

				throw __FILE__ ":build_db: Error reading Bloom filter CRC32";
			}
		}

		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] CRC32 data has been read." << endl;
		#endif // DEBUG_DB

		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read the Bloom filter info and number of filter bits
		vector<FilterInfo> filter_info(num_filter);

		for(size_t i = 0;i < num_filter;++i){

			binary_read(*fin_ptr[i], filter_info[i]);

			if( !(*fin_ptr[i]) ){

				#ifdef DEBUG_DB
				cerr << "[" << mpi_rank << "] Error reading Bloom filter info from " << m_bloom_files[i] << endl;
				#endif // DEBUG_DB

				throw __FILE__ ":build_db: Error reading Bloom filter info";
			}
		}

		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] All meta-data information has been read." << endl;
		#endif // DEBUG_DB

		// Before we start reading the Bloom filter bits and write the transposted data, we will open
		// the output file and write all of the header information.
		DBFileHeader output_header;
		
		output_header.crc32 = 0; // <-- This will be overwritten after we have transposed the Bloom filters
		output_header.kmer_len = m_param.kmer_len;
		output_header.num_hash = m_param.num_hash;
		output_header.log_2_filter_len = m_param.log_2_filter_len;
		output_header.num_filter = num_filter;

		// Not using compression for now (i.e. RLE_COMPRESSION). Compression is currently
		// too computationally intensive for SRA-scale database construction.
		output_header.compression = NO_COMPRESSION;

		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] About to open " << m_filename << " for database output" << endl;
		#endif // DEBUG_DB

		ofstream fout(m_filename.c_str(), ios::binary);

		if(!fout){
			throw __FILE__ ":build_db: Unable to open output file for writing";
		}

		// Write the header information for this file as a placeholder. We will need
		// to rewrite the header when we have calculated the crc32 value for all of 
		// the bitslices and the start of the metadata block.
		binary_write(fout, output_header);

		if(!fout){
			throw __FILE__ ":build_db: Error writing database file header placeholder";
		}

		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] Wrote output file header information." << endl;
		#endif // DEBUG_DB

		// Now we need to transpose the Bloom filter data (where filters {a,b,c,...} are
		// the source of the data):
		// |a|b|c|d|e|f|... <-- slice 0
		// |1|0|0|1|1|0|... <-- slice 1
		// |0|1|0|0|0|0|... <-- slice 2 
		// |1|0|1|1|0|1|... <-- slice 3
		// |0|1|1|0|0|0|... <-- slice 4
		// ...

		// The destination buffer of bitslices. To avoid using too much memory, we will write
		// the data in chunks of max_buffer_slice slices, each slice of length num_filter.
		//
		// Please note that max_buffer_slice must always be a multiple of BITS_PER_BLOCK (to ensure that
		// all binary reads are BLOCK aligned).
		// ** Lustre  (AWS FSX) has trouble with max_buffer_slice values that are too low **
		
		//const size_t max_buffer_slice = 64*BitVector::BITS_PER_BLOCK; <-- too small
		//const size_t max_buffer_slice = 512*BitVector::BITS_PER_BLOCK; <-- too small
		//const size_t max_buffer_slice = 262144*BitVector::BITS_PER_BLOCK; // 512 MB dest buffer
		const size_t max_buffer_slice = 524288*BitVector::BITS_PER_BLOCK; // 1024 MB dest buffer

		const size_t bytes_per_slice = num_filter/BitVector::BITS_PER_BLOCK + 
			( (num_filter%BitVector::BITS_PER_BLOCK) > 0 ? 1 : 0);

		unsigned char *dest = new unsigned char[max_buffer_slice*bytes_per_slice];

		if(dest == NULL){
			throw __FILE__ ":build_db: Unable to allocate destination buffer";
		}

		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] Transposing Bloom filters using a max_buffer_slice = " << max_buffer_slice << endl;
		#endif // DEBUG_DB

		// Start by reading num_buffer_slice bits at a time from each filter
		for(size_t i = 0;i < filter_len;i += max_buffer_slice){

			const size_t num_buffer_slice = min(max_buffer_slice, filter_len - i);

			// Unset all of the bits in the destination buffer so
			// we only need to set bits that are equal to one (rather than both setting *and* unsetting bits)
			const size_t curr_dest_len = num_buffer_slice*bytes_per_slice;

			memset(dest, 0, curr_dest_len);

			BitVector src(num_buffer_slice);

			// For each input bloom filter file:
			for(size_t j = 0;j < num_filter;++j){

				binary_read(*fin_ptr[j], src);

				if( !(*fin_ptr[j]) ){
					throw __FILE__ ":build_db: Error reading filter bytes";
				}

				// Keep a running update of the per-filter CRC32 values
				filter_crc32[j].second = ::crc32_z( filter_crc32[j].second, 
					src.ptr(), src.num_block() );
				
				// Copy the src slice from the j^th Bloom filter into the
				// j^th colum of destination buffer <-- this is the bitwise
				// transposition at the heart of the bitsliced Bloom filter
				// approach!
				BitVector tmp;

				for(size_t k = 0;k < num_buffer_slice;++k){

					if( src.get_bit(k) ){
						
						// The attach_bufer()/detach_buffer() functions allow the
						// BitVector set_bit() member function to write to an
						// arbitrary array of bytes
						tmp.attach_buffer(dest + bytes_per_slice*k, num_filter);

						tmp.set_bit(j);

						tmp.detach_buffer();
					}
				}
			}

			// Update the crc32 value with the contents of the current destination buffer
			output_header.crc32 = ::crc32_z( output_header.crc32, dest, curr_dest_len);

			// Write the destination buffer (i.e. the transpose of the Bloom filters) to disk
			fout.write( (const char*)dest, curr_dest_len );

			if(!fout){
				throw __FILE__ ":build_db: Unable to write transpose buffer to disk";
			}
		}

		if(dest != NULL){

			delete [] dest;
			dest = NULL;
		}

		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] Finished writing bitslices" << endl;
		#endif // DEBUG_DB

		// Manually deallocate (which will close) the input files. This is required by our use of pointers to ifstream
		// objects (which is, in turn, required to retain compatibility with older C++ compilers which do not have a
		// copy constructor defined for ifstream).
		for(vector<ifstream*>::iterator i = fin_ptr.begin();i != fin_ptr.end();++i){

			if(*i != NULL){

				delete *i;
				*i = NULL;
			}
		}

		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] Validating per-filter crc32 checksum values." << endl;
		#endif // DEBUG_DB

		// Validate the per-filter CRC32 values
		bool valid_filter_crc32 = true;

		for(size_t i = 0;i < num_filter;++i){

			if(filter_crc32[i].first != filter_crc32[i].second){

				valid_filter_crc32 = false;

				#ifdef DEBUG_DB
				cerr << "[" << mpi_rank << "] Error: Invalid CRC32 for Bloom filter " << m_bloom_files[i] << endl;
				cerr << "[" << mpi_rank << "] Read CRC32 = " << std::hex << filter_crc32[i].first 
					<< ", but computed CRC32 = " << filter_crc32[i].second << std::dec << endl;
				#endif // DEBUG_DB
			}
		}
		
		if(!valid_filter_crc32){
			throw __FILE__ ":build_db: One or more invalid Bloom filter CRC32 values";
		}

		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] Writing metadata." << endl;
		#endif // DEBUG_DB

		// We now need to write the filter metadata locations (i.e. the positions of each FilterInfo structure in the
		// file). Since we don't know the sizes yet, write a dummy array as a placeholder that we will
		// fill in later with the actual values.
		unsigned long int *info_loc_buffer = new unsigned long int [num_filter];
		
		if(info_loc_buffer == NULL){
			throw __FILE__ ":build_db: Unable to allocate metadata location buffer";
		}

		// Record the output file location of for the start of the metadata location array
		// (so we can fill this information in later)
		output_header.info_start = fout.tellp();
		
		memset( info_loc_buffer, 0, num_filter*sizeof(unsigned long int) );
			
		fout.write( (char*)info_loc_buffer, num_filter*sizeof(unsigned long int) );
			
		if(!fout){
			throw __FILE__ ":build_db: Error writing dummy metadata location buffer";
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// Write the Bloom filter info.
		for(size_t i = 0;i < num_filter;++i){

			// Record the starting file location of this filters metadata in the output file
			info_loc_buffer[i] = fout.tellp();

			binary_write(fout, filter_info[i]);

			if( !fout ){

				#ifdef DEBUG_DB
				cerr << "[" << mpi_rank << "] Error writing Bloom filter metadata from " << m_bloom_files[i] << endl;
				#endif // DEBUG_DB

				throw __FILE__ ":build_db: Error writing Bloom filter info";
			}
		}
		
		// Update the location array in the output file so we known where each filter metadata record starts
		fout.seekp(output_header.info_start);
		
		fout.write( (char*)info_loc_buffer, num_filter*sizeof(unsigned long int) );
			
		if(!fout){
			throw __FILE__ ":build_db: Error writing metadata location buffer";
		}
		
		delete [] info_loc_buffer;

		// Finally, we need to rewrite the header to update the crc32 value and the
		// start of the information block
		fout.seekp(0);

		binary_write(fout, output_header);

		if(!fout){
			throw __FILE__ ":build_db: Error writing database file header (final)";
		}

		fout.close();

		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] Finished building dabase file" << endl;
		#endif // DEBUG_DB
	}
	catch(const char *error){
		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] Caught the error " << error << endl;
		#endif // DEBUG_DB
		return false;
	}
	catch(const string error){
		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] Caught the error " << error << endl;
		#endif // DEBUG_DB
		return false;
	}
	catch(...){
		#ifdef DEBUG_DB
		cerr << "[" << mpi_rank << "] Caught an unhandled error" << endl;
		#endif // DEBUG_DB

		return false;
	}
	
	return true;
}
