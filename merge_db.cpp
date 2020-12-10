// Merge CALDERA database files to avoid partially filled files
// J. D. Gans
// Bioscience Division, B-11
// Los Alamos National Laboratory
// October 6, 2020

#include <iostream>
#include <iomanip>
#include <fstream>

#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <zlib.h>

#include "caldera.h"
#include "bloom.h"
#include "options.h"
#include "string_conversion.h"
#include "keys.h"
#include "file_util.h"

using namespace std;

pair<size_t , string> merge_database_files(const string &m_file_1, const string &m_file_2, const size_t &m_max_num_filters);

int main(int argc, char *argv[])
{
	try{
		
		// Command line arguments
		// <database file 1> <database file 2> ...
		// [-?|-h] (display command line options)

		deque<string> input_filename;
		
		const char* options = "h?";
		//int config_opt = 0;
		int long_index = 0;

		struct option long_opts[] = {
			{0,0,0,0} // Terminate options list
		};

		int opt_code;
		opterr = 0;

		bool print_usage = (argc == 1);
	
		while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

			switch( opt_code ){
				case 0:
					cerr << "Unknown flag!" << endl;
					break;
				case 'h':
				case '?':
					print_usage = true;
					break;
				default:
					cerr << '\"' << (char)opt_code << "\" is not a valid option!" << endl;
					break;
			};
		}

		for(int i = optind;i < argc;i++){
			input_filename.push_back(argv[i]);
		}

		if(print_usage){

			cerr << "Usage: "  << endl;
			cerr << "\t" << argv[0] << " <database file 1> <database file 2> ..." << endl;
			return EXIT_SUCCESS;
		}
		
		if(input_filename.size() < 2){

			cerr << "Please specify 2 or more database files to merge" << endl;
			return EXIT_SUCCESS;
		}

		// To keep the database filesizes reasonable, the maximum number of filters per database file
		// is adjusted based on the Bloom filter size. This is currently hard-coded, but could become
		// a user parameter

		MAP<size_t /*log 2 bloom filter size*/, size_t /* filters per database file*/> num_filters_per_database_file;

		for(size_t log_len = DEFAULT_MIN_LOG_2_FILTER_LEN;log_len <= DEFAULT_MAX_LOG_2_FILTER_LEN;++log_len){

			// If possible, we would like to have MAX_NUM_FILTER_CHUNK Bloom filters per database file (but no more).
			// However, we are not allowed to exceed MAX_DATABASE_FILE_SIZE_IN_GB gigabytes per database file!
			// In the equation for num bloom, 
			const size_t bits_per_byte = 8;
			const size_t num_bloom = (MAX_DATABASE_FILE_SIZE_IN_GB * bits_per_byte * GB)/(1UL << log_len);

			num_filters_per_database_file[log_len] = min(MAX_NUM_FILTER_CHUNK, num_bloom);
		}
		
		// There are several conditions that must be met for merging a pair of database files:
		// 		1) We will only attempt to merge database files that have the same Bloom filter parameters:
		//			kmer_len
		//			hash_func
		// 			num_hash
		//			filter_len
		//		2) Both files must have *less* than the maximum allowed number of Bloom filters

		MAP<string /*filename*/, DBFileHeader> header;

		// Group all of the input files by {kmer_len, hash_func, num_hash, filter_len}, as we can only merge databases
		// that share these key parameters.
		MULTIMAP<BloomParam, string /*filename*/> bloom_groups;

		for(deque<string>::const_iterator f = input_filename.begin();f != input_filename.end();++f){
			
			ifstream fin(f->c_str(), ios::binary);
			
			if(!fin){
			
				cerr << "Unable to open " << *f << " for reading" << endl;
				return EXIT_FAILURE;
			}
			
			DBFileHeader local_header;
			
			binary_read(fin, local_header);
			
			if(!fin){
			
				cerr << "Unable to read database header" << endl;
				return EXIT_FAILURE;
			}

			// Does the current database have any room?
			MAP<size_t /*log 2 bloom filter size*/, size_t /* filters per database file*/>::const_iterator max_iter = 
				num_filters_per_database_file.find(local_header.log_2_filter_len);

			if( max_iter == num_filters_per_database_file.end() ){

				cerr << "In database " << *f << ", could not lookup the maximum number of filters for a log2 filter length of " 
					<< local_header.log_2_filter_len << endl;

				return EXIT_FAILURE;
			}

			if(max_iter->second <= local_header.num_filter){

				//cerr << *f << " is full and not a merge candidate" << endl;
				continue;
			}

			// If we get here, then this database is a merge candidate. Make sure we don't have any duplicate
			// filenames
			if( header.find(*f) != header.end() ){
				
				cerr << *f << " appears more than once in the input file list" << endl;
				return EXIT_FAILURE;
			}

			header[*f] = local_header;

			BloomParam local_param;

			local_param.kmer_len = local_header.kmer_len;
			local_param.log_2_filter_len = local_header.log_2_filter_len;
			local_param.num_hash = local_header.num_hash;
			local_param.hash_func = local_header.hash_func;

			bloom_groups.insert( make_pair(local_param, *f) );
		}

		const vector<BloomParam> bloom_keys = keys(bloom_groups);
		const size_t num_bloom_groups = bloom_keys.size();

		cerr << "Found " << num_bloom_groups << " distinct Bloom parameter groups" << endl;

		for(size_t i = 0;i < num_bloom_groups;++i){

			typedef MULTIMAP<BloomParam, string /*filename*/>::const_iterator I;
			const pair<I, I> range = bloom_groups.equal_range(bloom_keys[i]);

			deque< pair<size_t /*number of filters*/, string /*filename*/> > db_files;

			for(I j = range.first;j != range.second;++j){

				// Look up the number of Bloom filters in this database file
				MAP<string /*filename*/, DBFileHeader>::const_iterator header_iter = header.find(j->second);

				if( header_iter == header.end() ){
					throw __FILE__ ":main: Unable to lookup header information";
				}

				db_files.push_back( make_pair(header_iter->second.num_filter, j->second) );
			}
			
			// Sort in ascending order by the number of filters
			sort( db_files.begin(), db_files.end() );

			cerr << "Bloom parameters for group " << i << " of " << num_bloom_groups << endl;
			cerr << "log_2_filter_len = " << bloom_keys[i].log_2_filter_len << endl;
			cerr << "num_hash = " << bloom_keys[i].num_hash << endl;

			MAP<size_t /*log 2 bloom filter size*/, size_t /* filters per database file*/>::const_iterator max_iter = 
				num_filters_per_database_file.find(bloom_keys[i].log_2_filter_len);

			if( max_iter == num_filters_per_database_file.end() ){

				cerr << "Could not lookup the maximum number of filters for a log2 filter length of " 
					<< bloom_keys[i].log_2_filter_len << endl;

				return EXIT_FAILURE;
			}

			// Keep merging pairs of files until we can't merge any more
			while(db_files.size() > 1){

				// The choice of which file is labeled file_large vs file_small only impacts the
				// case where there are too many Bloom filters to fit in a single file. In this case,
				// we want to keep the ordering of the larger file intact, on the off chance that these
				// Bloom filters are for related SRA records and therefor more compressible.

				// The files are sorted in ascending order, call the smaller file "file_small"
				const string file_small = db_files.front().second;

				db_files.pop_front();

				// The files are sorted in ascending order, call the larger file "file_large"
				const string file_large = db_files.front().second;
				
				db_files.pop_front();

				cerr << "\tmerging:" << '\n' 
					<< "\t\t" << file_small << '\n'
					<< "\t\t" << file_large << endl;

				// Merge the database files, remove the old files and write the new file(s) to file_large
				// and, if needed, file_small. If file_large becomes full and there are filters that still need to
				// be written, write them to file_small (and return the number of filters in file_small).
				// If file_large is *not* full, then return the number of filters in file_large
				const pair<size_t , string> remainder = merge_database_files(file_large, file_small, 
					max_iter->second /*max num filters*/);

				if(remainder.first > 0){

					db_files.push_back(remainder);
					sort( db_files.begin(), db_files.end() );
				}
			}
		}
	}
	catch(const char *error){
		cerr << "Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(const string error){
		cerr << "Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(...){
		cerr << "Caught an unhandled error" << endl;
		return EXIT_FAILURE;
	}
	
	
	return EXIT_SUCCESS;
}

pair<size_t , string> merge_database_files(const string &m_file_1, const string &m_file_2, const size_t &m_max_num_filters)
{
	const string tmp_suffix = ".tmp";

	pair<size_t , string> ret(size_t(0), "");

	ifstream fin_1(m_file_1.c_str(), ios::binary);
			
	if(!fin_1){
		throw __FILE__ ":merge_database_files: Unable to open database file 1";
	}
	
	ifstream fin_2(m_file_2.c_str(), ios::binary);
			
	if(!fin_2){
		throw __FILE__ ":merge_database_files: Unable to open database file 2";
	}

	DBFileHeader src_header_1;
	DBFileHeader src_header_2;
	
	binary_read(fin_1, src_header_1);
	
	if(!fin_1){
		throw __FILE__ ":merge_database_files: Error reading header 1";
	}
	
	binary_read(fin_2, src_header_2);
	
	if(!fin_2){
		throw __FILE__ ":merge_database_files: Error reading header 2";
	}

	cerr << "\t\t\tSrc 1 has " << src_header_1.num_filter << " Bloom filters" << endl;
	cerr << "\t\t\tSrc 2 has " << src_header_2.num_filter << " Bloom filters" << endl;

	cerr << "\t\t\tMax Bloom filters/file =  " << m_max_num_filters << endl;

	// Make sure that these file headers are indeed compatible
	if( (src_header_1.log_2_filter_len != src_header_2.log_2_filter_len) || 
		(src_header_1.num_hash != src_header_2.num_hash) ||
		(src_header_1.kmer_len != src_header_2.kmer_len) ||
		(src_header_1.hash_func != src_header_2.hash_func) ){
		
		throw __FILE__ ":merge_database_files: Incompatible database files";
	}

	// For now, compression is not allowed
	if( (src_header_1.compression != NO_COMPRESSION) || (src_header_2.compression != NO_COMPRESSION) ){
		throw __FILE__ ":merge_database_files: Compressed database files are not currently supported";
	}

	if(src_header_1.num_filter >= m_max_num_filters){
		throw __FILE__ ":merge_database_files: Database file 1 has more than expected filters";
	}

	if(src_header_2.num_filter >= m_max_num_filters){
		throw __FILE__ ":merge_database_files: Database file 2 has more than expected filters";
	}

	// If we have more filters than can fit into a single file, then we will need to write out
	// the "remainder" in a second file
	const bool has_remainder = (src_header_1.num_filter + src_header_2.num_filter) > m_max_num_filters;

	const string dst_file_1 = m_file_1 + tmp_suffix;
	const string dst_file_2 = has_remainder ? m_file_2 + tmp_suffix : "";

	ofstream fout_1;
	ofstream fout_2;

	// Make sure the temp files do not already exist
	if( is_file(dst_file_1) ){
		throw __FILE__ ":merge_database_files: Temp database file 1 already exists";
	}

	if( has_remainder && is_file(dst_file_2) ){
		throw __FILE__ ":merge_database_files: Temp database file 2 already exists";
	}

	// Open the output file(s)
	fout_1.open(dst_file_1.c_str(), ios::binary);

	if(!fout_1){
		throw __FILE__ ":merge_database_files: Unable to open new_file_1";
	}

	if(has_remainder){

		fout_2.open(dst_file_2.c_str(), ios::binary);

		if(!fout_2){
			throw __FILE__ ":merge_database_files: Unable to open new_file_2";
		}
	}

	DBFileHeader dst_header_1 = src_header_1;
	DBFileHeader dst_header_2 = src_header_2;

	// Reset the CRC32 and info start values of the new files. These will be overwritten at the end
	dst_header_1.crc32 = crc32_z(0L, Z_NULL, 0);
	dst_header_2.crc32 = crc32_z(0L, Z_NULL, 0);

	dst_header_1.info_start = 0;
	dst_header_2.info_start = 0;

	if(has_remainder){

		dst_header_1.num_filter = m_max_num_filters;
		dst_header_2.num_filter = (src_header_1.num_filter + src_header_2.num_filter) - m_max_num_filters;

		// We will overwrite the existing m_file_2 with the remainder
		ret = make_pair(dst_header_2.num_filter, m_file_2);
	}
	else{

		dst_header_1.num_filter = src_header_1.num_filter + src_header_2.num_filter;
		dst_header_2.num_filter = 0;

		// If the file that contains the merged Bloom filters is *not* full, then we need
		// to return it so that it can potentially be merged with another file.
		if(dst_header_1.num_filter < m_max_num_filters){
			
			// We will overwrite the existing m_file_1 with the merged Bloom filters
			ret = make_pair(dst_header_1.num_filter, m_file_1);
		}
	}

	// Write the header information for the output file(s) as a placeholder. We will need
	// to rewrite the headers when we have calculated the crc32 value for all of 
	// the bitslices and the start of the metadata block.
	binary_write(fout_1, dst_header_1);

	if(!fout_1){
		throw __FILE__ ":merge_database_files: Error writing database file 1 header placeholder";
	}

	if(has_remainder){

		binary_write(fout_2, dst_header_2);

		if(!fout_2){
			throw __FILE__ ":merge_database_files: Error writing database file 1 header placeholder";
		}
	}

	cerr << "\t\t\tDst 1 has " << dst_header_1.num_filter << " Bloom filters" << endl;

	if(has_remainder){
		cerr << "\t\t\tDst 2 has " << dst_header_2.num_filter << " Bloom filters" << endl;
	}

	const size_t num_bitslice = src_header_1.filter_len();

	// Prepare the input and output buffers:
	//  <---------1------->   <--2-->
	//	AAAAAAAAAAAAAAAAAAA + aaaaaaa [bitslice 0]
	//  BBBBBBBBBBBBBBBBBBB + bbbbbbb [bitslice 1]
	//  CCCCCCCCCCCCCCCCCCC + ccccccc [bitslice 2]
	//
	// Both files have the same number of bitslices (i.e. num_bitslice)

	const size_t bytes_per_slice_src_1 = src_header_1.num_filter/BitVector::BITS_PER_BLOCK + 
		( (src_header_1.num_filter%BitVector::BITS_PER_BLOCK) > 0 ? 1 : 0);

	const size_t bytes_per_slice_src_2 = src_header_2.num_filter/BitVector::BITS_PER_BLOCK + 
		( (src_header_2.num_filter%BitVector::BITS_PER_BLOCK) > 0 ? 1 : 0);

	const size_t bytes_per_slice_dst_1 = dst_header_1.num_filter/BitVector::BITS_PER_BLOCK + 
		( (dst_header_1.num_filter%BitVector::BITS_PER_BLOCK) > 0 ? 1 : 0);

	const size_t bytes_per_slice_dst_2 = dst_header_2.num_filter/BitVector::BITS_PER_BLOCK + 
		( (dst_header_2.num_filter%BitVector::BITS_PER_BLOCK) > 0 ? 1 : 0);

	// The reading and writing buffers are all sized to accomodate max_num_slice_per_buffer
	// slices
	const size_t max_num_slice_per_buffer = 1024;

	unsigned char* src_buffer_1 = new unsigned char[max_num_slice_per_buffer * bytes_per_slice_src_1];
	
	if(src_buffer_1 == NULL){
		throw __FILE__ ":merge_database_files: Unable to allocate source buffer 1";
	}

	unsigned char* src_buffer_2 = new unsigned char[max_num_slice_per_buffer * bytes_per_slice_src_2];
	
	if(src_buffer_2 == NULL){
		throw __FILE__ ":merge_database_files: Unable to allocate source buffer 2";
	}

	unsigned char* dst_buffer_1 = new unsigned char[max_num_slice_per_buffer * bytes_per_slice_dst_1];
	
	if(dst_buffer_1 == NULL){
		throw __FILE__ ":merge_database_files: Unable to allocate destination buffer 1";
	}

	unsigned char* dst_buffer_2 = NULL;
	
	if(has_remainder){

		dst_buffer_2 = new unsigned char[max_num_slice_per_buffer * bytes_per_slice_dst_2];
	
		if(dst_buffer_2 == NULL){
			throw __FILE__ ":merge_database_files: Unable to allocate destination buffer 2";
		}
	}

	// Initialize the CRC32 values for the source database files (to verify bitsice integrity)
	uint32_t crc32_src_1 = crc32_z(0L, Z_NULL, 0);
	uint32_t crc32_src_2 = crc32_z(0L, Z_NULL, 0);

	// How many bits from each slice in src file 2 will be merged into dst file 1?
	const size_t num_src_2_bits_merge_1 = src_header_2.num_filter - dst_header_2.num_filter;

	// DEBUG
	//cerr << "About to merge bitslices; num_src_2_bits_merge_1 = " << num_src_2_bits_merge_1 << endl;
	//cerr << (has_remainder ? "Has remainder" : "Does *not* have remainder") << endl;
	//cerr << "num_bitslice = " << num_bitslice << endl;
	//cerr << "bytes_per_slice_src_1 = " << bytes_per_slice_src_1 << endl;
	//cerr << "bytes_per_slice_src_2 = " << bytes_per_slice_src_2 << endl;

	// Read and write the database files in chunks of max_num_slice_per_buffer bitslices
	for(size_t i = 0;i < num_bitslice;i += max_num_slice_per_buffer){

		const size_t num_buffer_slice = ( (i + max_num_slice_per_buffer) > num_bitslice ) ?
			 (num_bitslice - i) : max_num_slice_per_buffer;

		const size_t num_src_1_bytes = num_buffer_slice*bytes_per_slice_src_1;
		const size_t num_src_2_bytes = num_buffer_slice*bytes_per_slice_src_2;
		
		const size_t num_dst_1_bytes = num_buffer_slice*bytes_per_slice_dst_1;
		const size_t num_dst_2_bytes = num_buffer_slice*bytes_per_slice_dst_2;

		// Read from the source files into the source buffers
		fin_1.read( (char*)src_buffer_1, num_src_1_bytes );

		if(!fin_1){
			throw __FILE__ ":merge_database_files: Error reading bitslices from source file 1";
		}

		fin_2.read( (char*)src_buffer_2, num_src_2_bytes );

		if(!fin_2){
			throw __FILE__ ":merge_database_files: Error reading bitslices from source file 2";
		}

		// Update the CRC32 values as we go
		crc32_src_1 = ::crc32_z(crc32_src_1, src_buffer_1, num_src_1_bytes);
		crc32_src_2 = ::crc32_z(crc32_src_2, src_buffer_2, num_src_2_bytes);

		// Zero the destination buffers to make merging simple
		memset(dst_buffer_1, 0, num_dst_1_bytes);

		if(has_remainder){
			memset(dst_buffer_2, 0, num_dst_2_bytes);
		}

		// Merge the source buffers into the destination buffers
		for(size_t j = 0;j < num_buffer_slice;++j){

			// Copy the src 1 slice into dst 1
			memcpy(dst_buffer_1 + j*bytes_per_slice_dst_1, src_buffer_1 + j*bytes_per_slice_src_1, 
				bytes_per_slice_src_1);

			// Merging the src 2 slice is a little tricker, as we need to operate on individual bits.
			BitVector dst_1;
			BitVector src_2;

			dst_1.attach_buffer(dst_buffer_1 + j*bytes_per_slice_dst_1, dst_header_1.num_filter);
			src_2.attach_buffer(src_buffer_2 + j*bytes_per_slice_src_2, src_header_2.num_filter);

			for(size_t k = 0;k < num_src_2_bits_merge_1;++k){

				if( src_2.get_bit(k) ){
					dst_1.set_bit(src_header_1.num_filter + k);
				}
			}

			if(has_remainder){

				BitVector dst_2;

				dst_2.attach_buffer(dst_buffer_2 + j*bytes_per_slice_dst_2, dst_header_2.num_filter);

				for(size_t k = 0;k < dst_header_2.num_filter;++k){

					if( src_2.get_bit(num_src_2_bits_merge_1 + k) ){
						dst_2.set_bit(k);
					}
				}
				
				dst_2.detach_buffer();
			}

			dst_1.detach_buffer();
			src_2.detach_buffer();
		}

		// Write the destination buffers
		fout_1.write( (char*)dst_buffer_1, num_dst_1_bytes );

		if(!fout_1){
			throw __FILE__ ":merge_database_files: Error writing bitslices to destination file 1";
		}

		dst_header_1.crc32 = ::crc32_z(dst_header_1.crc32, dst_buffer_1, num_dst_1_bytes);

		if(has_remainder){

			fout_2.write( (char*)dst_buffer_2, num_dst_2_bytes );

			if(!fout_2){
				throw __FILE__ ":merge_database_files: Error writing bitslices to destination file 2";
			}
			
			dst_header_2.crc32 = ::crc32_z(dst_header_2.crc32, dst_buffer_2, num_dst_2_bytes);
		}
	}

	// DEBUG
	//cerr << "About to free buffers" << endl;

	// Clear the bitslice buffers
	if(src_buffer_1 != NULL){

		delete [] src_buffer_1;
		src_buffer_1 = NULL;
	}

	if(src_buffer_2 != NULL){

		delete [] src_buffer_2;
		src_buffer_2 = NULL;
	}

	if(dst_buffer_1 != NULL){

		delete [] dst_buffer_1;
		dst_buffer_1 = NULL;
	}

	if(crc32_src_1 != src_header_1.crc32){
		throw __FILE__ ":merge_database_files: Invalid CRC32 value for source database file 1";
	}

	if(crc32_src_2 != src_header_2.crc32){
		throw __FILE__ ":merge_database_files: Invalid CRC32 value for source database file 2";
	}

	// DEBUG
	//cerr << "About to write info placeholders" << endl;

	// Record the start of the information section
	dst_header_1.info_start = fout_1.tellp();

	if(has_remainder){
		dst_header_2.info_start = fout_2.tellp();
	}

	// We now need to write the filter metadata locations (i.e. the positions of each FilterInfo structure in the
	// file). Since we don't know the sizes yet, write a dummy array as a placeholder that we will
	// fill in later with the actual values.
	unsigned long int *info_loc_buffer_1 = new unsigned long int [dst_header_1.num_filter];
	
	if(info_loc_buffer_1 == NULL){
		throw __FILE__ ":merge_database_files: Unable to allocate metadata location buffer for destination database file 1";
	}
	
	memset( info_loc_buffer_1, 0, dst_header_1.num_filter*sizeof(unsigned long int) );
			
	fout_1.write( (char*)info_loc_buffer_1, dst_header_1.num_filter*sizeof(unsigned long int) );
		
	if(!fout_1){
		throw __FILE__ ":merge_database_files: Error writing dummy metadata location buffer to destination databaes file 1";
	}
	
	unsigned long int *info_loc_buffer_2 = NULL;
	
	if(has_remainder){

		info_loc_buffer_2 = new unsigned long int [dst_header_2.num_filter];
		
		if(info_loc_buffer_2 == NULL){
			throw __FILE__ ":merge_database_files: Unable to allocate metadata location buffer for destination database file 2";
		}
		
		memset( info_loc_buffer_2, 0, dst_header_2.num_filter*sizeof(unsigned long int) );
				
		fout_2.write( (char*)info_loc_buffer_2, dst_header_2.num_filter*sizeof(unsigned long int) );
			
		if(!fout_2){
			throw __FILE__ ":merge_database_files: Error writing dummy metadata location buffer to destination databaes file 2";
		}
	}

	// DEBUG
	//cerr << "About to jump to input info sections" << endl;

	// Skip the info_loc_buffer sections in the input files. We don't need this, as we are reading 
	// the information sequentially.
	fin_1.seekg( size_t( fin_1.tellg() ) + sizeof(unsigned long int)*src_header_1.num_filter);
	fin_2.seekg( size_t( fin_2.tellg() ) + sizeof(unsigned long int)*src_header_2.num_filter);

	// DEBUG
	//cerr << "About to write file 1 info section" << endl;

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// Write the Bloom filter info from the first source file
	for(size_t i = 0;i < src_header_1.num_filter;++i){

		FilterInfo info;

		binary_read(fin_1, info);

		// Record the starting file location of this filters metadata in the output file
		info_loc_buffer_1[i] = fout_1.tellp();

		binary_write(fout_1, info);

		if( !fout_1 ){
			throw __FILE__ ":merge_database_files: Error writing Bloom filter info to destination file 1";
		}
	}
	
	// DEBUG
	//cerr << "About to write file 2 info section" << endl;

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// Write the Bloom filter info from the second source file
	for(size_t i = 0;i < num_src_2_bits_merge_1;++i){

		FilterInfo info;

		binary_read(fin_2, info);

		// Record the starting file location of this filters metadata in the output file
		info_loc_buffer_1[src_header_1.num_filter + i] = fout_1.tellp();

		binary_write(fout_1, info);

		if( !fout_1 ){
			throw __FILE__ ":merge_database_files: Error writing Bloom filter info to destination file 1 (b)";
		}
	}

	if(has_remainder){

		for(size_t i = 0;i < dst_header_2.num_filter;++i){

			FilterInfo info;

			binary_read(fin_2, info);

			// Record the starting file location of this filters metadata in the output file
			info_loc_buffer_2[i] = fout_2.tellp();

			binary_write(fout_2, info);

			if( !fout_2 ){
				throw __FILE__ ":merge_database_files: Error writing Bloom filter info to destination file 2";
			}
		}
	}

	// DEBUG
	//cerr << "About to write file 1 header" << endl;

	// Rewind to write the header and info_loc information for both destination files
	fout_1.seekp(0);

	binary_write(fout_1, dst_header_1);

	if(!fout_1){
		throw __FILE__ ":merge_database_files: Error writing database file 1 header (final)";
	}

	// DEBUG
	//cerr << "About to write file 1 info start" << endl;

	fout_1.seekp(dst_header_1.info_start);
	
	fout_1.write( (char*)info_loc_buffer_1, dst_header_1.num_filter*sizeof(unsigned long int) );
		
	if(!fout_1){
		throw __FILE__ ":merge_database_files: Error writing metadata location buffer to destination file 1";
	}

	fout_1.close();

	// DEBUG
	//cerr << "Closed file 1" << endl;

	if(has_remainder){

		fout_2.seekp(0);

		binary_write(fout_2, dst_header_2);

		if(!fout_2){
			throw __FILE__ ":merge_database_files: Error writing database file 2 header (final)";
		}
		
		fout_2.seekp(dst_header_2.info_start);
	
		fout_2.write( (char*)info_loc_buffer_2, dst_header_2.num_filter*sizeof(unsigned long int) );
			
		if(!fout_2){
			throw __FILE__ ":merge_database_files: Error writing metadata location buffer to destination file 2";
		}

		fout_2.close();
	}

	if(info_loc_buffer_1 != NULL){

		delete [] info_loc_buffer_1;
		info_loc_buffer_1 = NULL;
	}

	if(info_loc_buffer_2 != NULL){

		delete [] info_loc_buffer_2;
		info_loc_buffer_2 = NULL;
	}

	// Overwrite the input files with the temporary files. Note that rename() only works for files
	// in the same filesystem.
	if(rename( dst_file_1.c_str(), m_file_1.c_str() ) != 0){

		cerr << "Unable to move first temp file " << dst_file_1 << " to " << m_file_1 << endl;
		throw __FILE__ ":merge_database_files: Error renaming database file 1";
	}

	if(has_remainder){

		if(rename( dst_file_2.c_str(), m_file_2.c_str() ) != 0){

			cerr << "Unable to move second temp file " << dst_file_2 << " to " << m_file_2 << endl;
			throw __FILE__ ":merge_database_files: Error renaming database file 2";
		}
	}
	else{

		// If there is *no* remainder, then we have merged all of the Bloom filters from m_file_2 into m_file_1 and
		// it is safe to delete m_file_2
		if(unlink( m_file_2.c_str() ) != 0){

			cerr << "Unable to remove file " << m_file_2 << endl;
			throw __FILE__ ":merge_database_files: Error removing database file 2";
		}
	}

	return ret;
}
