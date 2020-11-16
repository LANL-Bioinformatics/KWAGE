// Compress BIGSI++ database files
// J. D. Gans
// Bioscience Division, B-11
// Los Alamos National Laboratory
// October 26, 2020

#include <iostream>
#include <iomanip>
#include <fstream>
#include <unordered_set>
#include <queue>

#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <zlib.h>

#include "bigsi++.h"
#include "bloom.h"
#include "options.h"
#include "string_conversion.h"
#include "keys.h"
#include "file_util.h"
#include "slice_z.h"

using namespace std;

typedef size_t PayloadBuffer;
typedef unsigned int CountBuffer;

void compress_database_file(const string &m_file);
float bloom_distance(const BitVector &m_a, const BitVector &m_b);
unsigned int compress_rle(unsigned char** m_out, const unsigned char* m_in, 
	const unsigned int &m_in_bits,
	const unsigned int &m_num_payload_bits, 
	const unsigned int &m_num_count_bits, 
	const unsigned int &m_min_count);

void decode_rle(unsigned char* m_out, const unsigned int &m_out_bits, 
	const unsigned char* m_in, 
	const unsigned int &m_num_payload_bits, 
	const unsigned int &m_num_count_bits, 
	const unsigned int &m_min_count);

// Compute the size needed to store the compressed output, but do not acutally compress
unsigned int compress_rle(const unsigned char* m_in, const unsigned int &m_len,
	const unsigned int &m_num_payload_bits, const unsigned int &m_num_count_bits, const unsigned int &m_min_count);

int main(int argc, char *argv[])
{
	try{
		
		#ifdef NOT_NOW
		const unsigned int debug_num_payload_bits = 8;
		const unsigned int debug_num_count_bits = 4;
		const unsigned int debug_min_count = 4;

		for(unsigned int len = 5;len < 70;++len){
		//for(unsigned int len = 6;len < 7;++len){

			for(unsigned int trial = 0;trial < 10000;++trial){

				//cerr << "*** " << len << " ***" << endl;

				unsigned char *input = new unsigned char[len];

				for(unsigned int i = 0;i < len;++i){
					input[i] = rand()%256;
					//input[i] = 0;
				}

				// DEBUG
				//input[0] = 8;
				//input[1] = 33;
				//input[2] = 36;
				//input[3] = 132;
				//input[4] = 173;
				//input[5] = 195;

				sort(input, input + len);

				//cerr << "Input =";
			
				//for(unsigned int i = 0;i < len;++i){
				//	cerr << ' ' << int(input[i]);
				//}

				//cerr << endl;

				//cerr << "input bytes = " << len << endl;
				//cerr << "input bits = " << len*8 << endl;

				unsigned char* output = NULL;

				const unsigned int compressed_bits = 
					compress_rle( &output, (const unsigned char*)input, 8*len,
					debug_num_payload_bits, debug_num_count_bits, debug_min_count);

				//cerr << "compressed bits/bytes = " << compressed_bits << " / " << 
				//	( compressed_bits/8 + (compressed_bits%8 > 0) ) << endl;

				//throw "DEBUG";

				if(compressed_bits > 0){

					cerr << "[" << len << "] compressed bits/bytes = " << compressed_bits << " / " << 
						( compressed_bits/8 + (compressed_bits%8 > 0) ) << endl;

					cout << len << '\t' << compressed_bits << '\t' << 
						compress_rle((const unsigned char*)input, len,
							debug_num_payload_bits, debug_num_count_bits, debug_min_count)
						<< endl;

					unsigned char test[len];

					decode_rle(test, len*8, output, debug_num_payload_bits, debug_num_count_bits, debug_min_count);

					bool valid = true;

					for(unsigned int i = 0;i < len;++i){
						if(test[i] != input[i]){
							valid = false;
						}
					}

					if(!valid){

						cerr << "Failure for len = " << len << endl;

						for(unsigned int i = 0;i < len;++i){
							cout << i << '\t' << int(input[i]) << '\t' << int(test[i]) << endl;
						}

						throw "Invalid";
					}
				}
				//else{
				//	cerr << "Unable to compress" << endl;
				//}

				if(output != NULL){

					delete [] output;
					output = NULL;
				}

				delete [] input;
			}
		}

		throw "DEBUG";
		#endif // NOT_NOW

		time_t profile = time(NULL);

		// Command line arguments
		// <database file 1> ...
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
			cerr << "\t" << argv[0] << " <database file 1> ..." << endl;
			return EXIT_SUCCESS;
		}
		
		if( input_filename.empty() ){

			cerr << "Please specify at least one database file to compress" << endl;
			return EXIT_SUCCESS;
		}
		
		for(deque<string>::const_iterator f = input_filename.begin();f != input_filename.end();++f){
			compress_database_file(*f);
		}

		profile = time(NULL) - profile;

		cerr << "Compressed in " << profile << " src" << endl;
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

void compress_database_file(const string &m_file)
{
	ifstream fin(m_file.c_str(), ios::binary);
			
	if(!fin){
		throw __FILE__ ":compress_database_file: Unable to open database file";
	}
	
	DBFileHeader src_header;
	
	binary_read(fin, src_header);
	
	if(!fin){
		throw __FILE__ ":compress_database_file: Error reading header";
	}

	cerr << "Attempting to compress " << m_file << endl;
	cerr << "\tInput DB has " << src_header.num_filter << " Bloom filters" << endl;

	// Make sure that these file headers are indeed compatible
	if( src_header.compression != NO_COMPRESSION){

		cerr << "\tDB file is already compressed" << endl;
		return;
	}

	const string dst_file = m_file + "z";

	ofstream fout;

	// Make sure the output file does not already exist
	if( is_file(dst_file) ){
		throw __FILE__ ":compress_database_file: Compressed file already exists";
	}

	// Open the output file(s)
	fout.open(dst_file.c_str(), ios::binary);

	if(!fout){
		throw __FILE__ ":compress_database_file: Unable to open compressed file";
	}

	DBFileHeader dst_header = src_header;

	dst_header.compression = RLE_COMPRESSION;

	// Reset the info start values of the new files. These will be overwritten at the end
	dst_header.info_start = 0;

	// Write the header information for the output file(s) as a placeholder. We will need
	// to rewrite the headers when we have calculated the crc32 value for all of 
	// the bitslices and the start of the metadata block.
	binary_write(fout, dst_header);

	if(!fout){
		throw __FILE__ ":compress_database_file: Error writing compressed database header placeholder";
	}

	const size_t num_bitslice = src_header.filter_len();

	cerr << "\tInput DB has " << num_bitslice << " bitslices (log_2 = " 
		<< src_header.log_2_filter_len << ")" << endl;

	// The RLE_COMPRESSION format requires the following header:
	//	1 byte for the num_payload_bits
	//	1 byte for the num_count_bits
	//	1 byte for the min_count
	//	num_bitslice bytes for the length of each compressed buffer

	unsigned char *compressed_bytes = new unsigned char [num_bitslice];

	if(compressed_bytes == NULL){
		throw __FILE__ ":compress_database_file: Unable to allocate compressed bytes";
	}

	// Even though we don't yet know the number of bytes in each compressed slice, write this array
	// as a placeholder.
	fout.write( (char*)compressed_bytes, num_bitslice);

	const size_t bytes_per_slice_src = src_header.num_filter/BitVector::BITS_PER_BLOCK + 
		( (src_header.num_filter%BitVector::BITS_PER_BLOCK) > 0 ? 1 : 0);

	// The reading and writing buffers are all sized to accomodate max_num_slice_per_buffer
	// slices
	const size_t max_num_slice_per_buffer = 2048;
	
	// The size of both the source and destination buffers
	const size_t buffer_size = max_num_slice_per_buffer * bytes_per_slice_src;

	unsigned char* src_buffer = new unsigned char[buffer_size];
	
	if(src_buffer == NULL){
		throw __FILE__ ":compress_database_file: Unable to allocate source buffer";
	}

	unsigned char* dst_buffer = new unsigned char[buffer_size];
	
	if(dst_buffer == NULL){
		throw __FILE__ ":compress_database_file: Unable to allocate destination buffer";
	}

	// Initialize the CRC32 values for the source database file (to verify bitsice integrity)
	uint32_t crc32_src = crc32_z(0L, Z_NULL, 0);

	// TODO: How to consistently order the filters to enable maximum compression?
	// Scheme #1: Randomly sample a small number of bit slices and sort the filters by
	// 			  the average bit value (average over sampled slices).
	// Scheme #1: Randomly sample a small number of bit slices and permute order (Monte Carlo?,
	//			  greedy?, other?) to optimize the compression ratio.

	// Remember where the source data starts
	const size_t src_data_start = fin.tellg();

	fin.read( (char*)src_buffer, buffer_size );

	if(!fin){
		throw __FILE__ ":compress_database_file: Error reading bitslices from source file";
	}

	///////////////////////////////////////////////////////////////////////////////////////
	#define ORDER_BLOOM_FILTERS
	#ifdef ORDER_BLOOM_FILTERS
	
	cerr << "\tOrdering Bloom filters ... ";

	// In order to compute the pairwise distances, we need to transpose the bitslices into
	// Bloom filters
	vector< BitVector > bloom( src_header.num_filter, BitVector(max_num_slice_per_buffer) );

	// The average Bloom filter is the center of the cluster
	vector<float> ave_bloom(max_num_slice_per_buffer);

	for(size_t i = 0;i < src_header.num_filter;++i){
		bloom[i].unset_all_bits();
	}

	// Find the center of the partial Bloom filters
	for(size_t i = 0;i < max_num_slice_per_buffer;++i){

		BitVector src;

		src.attach_buffer(src_buffer + i*bytes_per_slice_src, src_header.num_filter);

		for(unsigned int j = 0;j < src_header.num_filter;++j){
			
			if( src.get_bit(j) ){

				bloom[j].set_bit(i);

				ave_bloom[i] += 1.0f;
			}
		}

		src.detach_buffer();
	}

	// Normalization the average partial Bloom filter
	for(size_t i = 0;i < max_num_slice_per_buffer;++i){
		ave_bloom[i] /= src_header.num_filter;
	}

	// Find the Bloom filter that is *farthest* from the center of the Bloom filters
	float max_dist = -1.0;
	size_t index = 0;

	for(size_t i = 0;i < src_header.num_filter;++i){

		float dist = 0.0;

		for(size_t j = 0;j < max_num_slice_per_buffer;++j){
			
			if( bloom[i].get_bit(j) ){

				const float delta = 1.0f - ave_bloom[j];

				dist += delta*delta;
			}
			else{
				dist += ave_bloom[j]*ave_bloom[j];
			}
		}

		if(dist > max_dist){

			max_dist = dist;
			index = i;
		}
	}

	if(max_dist < 0.0){
		throw __FILE__ ":compress_database_file: Unable to find max_index";
	}

	vector<bool> filter_pool(src_header.num_filter); // <-- which filters have already been ordered?
	vector<unsigned int> filter_order(src_header.num_filter);

	size_t curr_col = 0;

	// The max_index Bloom filter is the first column
	filter_order[index] = curr_col;
	++curr_col;

	filter_pool[index] = true;

	// Iterate until we have ordered all Bloom filters
	for(;curr_col < src_header.num_filter;++curr_col){

		float min_dist = 1.0e9;
		size_t next_index = 0;

		for(size_t i = 0;i < src_header.num_filter;++i){

			// Skip filters that have already been ordered
			if(filter_pool[i]){
				continue;
			}

			const float dist = bloom_distance(bloom[i], bloom[index]);

			if(dist < min_dist){

				min_dist = dist;
				next_index = i;
			}
		}

		index = next_index;

		filter_order[index] = curr_col;
		filter_pool[index] = true;
	}
	
	cerr << "done." << endl;

	#else // Default order

	vector<unsigned int> filter_order(src_header.num_filter);

	for(size_t i = 0;i < src_header.num_filter;++i){
		filter_order[i] = i;
	}

	#endif // ORDER_BLOOM_FILTERS

	///////////////////////////////////////////////////////////////////////////////////////
	// Find the optimal RLE parameters
	unsigned int best_num_payload_bits = 0;
	unsigned int best_num_count_bits = 0;
	unsigned int best_min_count = 0;
	unsigned int best_rle_size = 0xFFFFFFFF;

	cerr << "\tOptimzing RLE parameters" << endl;

	for(unsigned int num_payload_bits = 1;num_payload_bits < 8;++num_payload_bits){

		for(unsigned int num_count_bits = 1;num_count_bits < 8;++num_count_bits){
			
			for(unsigned int min_count = 1;min_count < 8;++min_count){

				unsigned int compressed_bytes = 0;

				for(size_t i = 0;i < max_num_slice_per_buffer;++i){

					BitVector src;

					src.attach_buffer(src_buffer + i*bytes_per_slice_src, src_header.num_filter);

					BitVector dst(src_header.num_filter);

					dst.unset_all_bits();

					// Reorder the filter bits
					for(unsigned int k = 0;k < src_header.num_filter;++k){

						if( src.get_bit(k) ){
							dst.set_bit(filter_order[k]);
						}
					}

					src.detach_buffer();

					const unsigned int rle_bytes = compress_rle(dst.ptr(), bytes_per_slice_src,
						num_payload_bits, num_count_bits, min_count);

					compressed_bytes += min(rle_bytes, (unsigned int)bytes_per_slice_src);
				}

				if(best_rle_size > compressed_bytes){

					best_num_payload_bits = num_payload_bits;
					best_num_count_bits = num_count_bits;
					best_min_count = min_count;

					best_rle_size = compressed_bytes;
				}
			}
		}
	}

	cerr << "\t\tBest RLE compression = " << float(best_rle_size)/buffer_size << endl;
	cerr << "\t\tbest_num_payload_bits = " << best_num_payload_bits << endl;
	cerr << "\t\tbest_num_count_bits = " << best_num_count_bits << endl;
	cerr << "\t\tbest_min_count = " << best_min_count << endl;

	///////////////////////////////////////////////////////////////////////////////////////

	vector<size_t> hist(256);
	size_t total = 0;

	// Can we further compress the data with Huffman or arthimetic encoding?
	for(size_t i = 0;i < max_num_slice_per_buffer;++i){

		BitVector src;

		src.attach_buffer(src_buffer + i*bytes_per_slice_src, src_header.num_filter);

		BitVector dst(src_header.num_filter);

		dst.unset_all_bits();

		// Reorder the filter bits
		for(unsigned int k = 0;k < src_header.num_filter;++k){

			if( src.get_bit(k) ){
				dst.set_bit(filter_order[k]);
			}
		}

		src.detach_buffer();

		// DEBUG
		for(unsigned int j = 0;j < bytes_per_slice_src;++j){

			++hist[ dst.ptr()[j] ];
			++total;

			//++hist[ zbuffer[j] & 0xF];
			//++total;
			//++hist[ (zbuffer[j] >> 4)  & 0xF];
			//++
		}

		continue;

		unsigned char *zbuffer = NULL;

		unsigned int rle_bits = compress_rle(&zbuffer, 
			dst.ptr(), 8*bytes_per_slice_src /*bits not bytes*/,
			best_num_payload_bits, best_num_count_bits, best_min_count);

		if( ( rle_bits/8 + (rle_bits%8 > 0) ) == bytes_per_slice_src ){

			delete [] zbuffer;
			zbuffer = NULL;

			rle_bits = 0;
		}

		if(rle_bits > 0){

			unsigned int bytes = rle_bits/8 + (rle_bits%8 > 0);

			for(unsigned int j = 0;j < bytes;++j){

				++hist[ zbuffer[j] ];
				++total;

				//++hist[ zbuffer[j] & 0xF];
				//++total;
				//++hist[ (zbuffer[j] >> 4)  & 0xF];
				//++total;
			}

			delete [] zbuffer;
			zbuffer = NULL;
		}
		else{

			//const unsigned char* buffer = dst.ptr();

			//for(unsigned int j = 0;j < bytes_per_slice_src;++j){
			//	++hist[ buffer[j] ];
			//}
		}
	}

	priority_queue<float, std::vector<float>, std::greater<float> > tree;
	float shannon_entropy = 0.0;

	cout << "Byte histogram:" << endl;

	for(unsigned int i = 0;i < 256;++i){

		cout << i << '\t' << hist[i] << endl;

		const float w = float(hist[i])/total;

		tree.push( float(w) );

		if(w > 0.0){
			shannon_entropy += -w*log(w)/log(2.0);
		}
	}

	const float bit_norm = 8.0;

	cerr << "Shannon entropy = " << shannon_entropy << " (ideal = " << shannon_entropy/bit_norm << ")" << endl;

	float bits_per_symbol = 0.0;

	while(tree.size() > 1){

		float node = tree.top();
		bits_per_symbol += tree.top();

		tree.pop();

		node += tree.top();
		bits_per_symbol += tree.top();

		tree.pop();

		tree.push(node);
	}

	cerr << "Huffman coding entropy = " << bits_per_symbol << " (ideal = " << bits_per_symbol/bit_norm << ")" << endl;

	///////////////////////////////////////////////////////////////////////////////////////

	CompressSlice<MAX_COMPRESSED_BYTES> zlib_buffer;

	// Rewind the source data file
	fin.seekg(src_data_start);

	size_t uncompressed_summary = 0;
	size_t zlib_compressed_summary = 0;
	size_t rle_compressed_summary = 0;
	unsigned int min_rle_len = 0xFFFFFFFF;
	unsigned int max_rle_len = 0x0;

	// Read and write the database files in chunks of max_num_slice_per_buffer bitslices
	for(size_t i = 0;i < num_bitslice;i += max_num_slice_per_buffer){

		const size_t num_buffer_slice = ( (i + max_num_slice_per_buffer) > num_bitslice ) ?
			 (num_bitslice - i) : max_num_slice_per_buffer;

		const size_t num_src_bytes = num_buffer_slice*bytes_per_slice_src;
		
		uncompressed_summary += num_src_bytes;

		// Read from the source files into the source buffers
		fin.read( (char*)src_buffer, num_src_bytes );

		if(!fin){
			throw __FILE__ ":compress_database_file: Error reading bitslices from source file";
		}

		// Update the CRC32 values as we go
		crc32_src = ::crc32_z(crc32_src, src_buffer, num_src_bytes);

		unsigned char* dst_head = dst_buffer;
		size_t dst_size = 0;

		// Attempt to compress the source buffer and store the results in the destination buffer
		for(size_t j = 0;j < num_buffer_slice;++j){

			BitVector src;

			src.attach_buffer(src_buffer + j*bytes_per_slice_src, src_header.num_filter);

			BitVector dst(src_header.num_filter);

			dst.unset_all_bits();

			// Reorder the filter bits
			for(unsigned int k = 0;k < src_header.num_filter;++k){

				if( src.get_bit(k) ){
					dst.set_bit(filter_order[k]);
				}
			}

			src.detach_buffer();

			unsigned char* zbuffer = NULL;

			// bit-wise RLE compression
			unsigned int rle_bits = compress_rle(&zbuffer, 
				dst.ptr(), 8*bytes_per_slice_src /*bits, not bytes*/,
				best_num_payload_bits, best_num_count_bits, best_min_count);

			// zlib compression (for comparison only)
			if( zlib_buffer.compress(dst.ptr(), bytes_per_slice_src) ){
				zlib_compressed_summary += zlib_buffer.size();
			}
			else{
				zlib_compressed_summary += bytes_per_slice_src;
			}

			// Compression that only improves the result by a fraction of a byte
			// is not worth the effort! In addition, we can not store the value 256
			// in an unsigned char!
			if( ( rle_bits/8 + (rle_bits%8 > 0) ) == bytes_per_slice_src ){

				delete [] zbuffer;
				zbuffer = NULL;

				rle_bits = 0;
			}

			if(rle_bits > 0){

				compressed_bytes[i + j] = rle_bits/8 + (rle_bits%8 > 0);

				rle_compressed_summary += compressed_bytes[i + j];

				min_rle_len = min(min_rle_len, (unsigned int)compressed_bytes[i + j]);
				max_rle_len = max(max_rle_len, (unsigned int)compressed_bytes[i + j]);

				dst_size += compressed_bytes[i + j];

				if(dst_size > buffer_size){
					throw __FILE__ ":compress_database_file: Destination buffer overflow";
				}

				memcpy( dst_head, zbuffer, compressed_bytes[i + j] );

				dst_head += compressed_bytes[i + j];

				// We only need to free the RLE zbuffer when rle_bits > 0
				delete [] zbuffer;
				zbuffer = NULL;
			}
			else{ // Failed to compress

				rle_compressed_summary += bytes_per_slice_src;

				//min_rle_len = min(min_rle_len, (unsigned int)bytes_per_slice_src);
				//max_rle_len = max(max_rle_len, (unsigned int)bytes_per_slice_src);

				compressed_bytes[i + j] = 0; // <-- zero is a special value that means no compression

				dst_size += bytes_per_slice_src;

				if(dst_size > buffer_size){
					throw __FILE__ ":compress_database_file: Destination buffer overflow";
				}

				// Make sure to store the uncompress, but re-ordered bitslice
				memcpy( dst_head, dst.ptr(), bytes_per_slice_src );
				
				dst_head += bytes_per_slice_src;
			}

			// DEBUG
			//cout << i + j << '\t' << (int)compressed_bytes[i + j] << " / " << (int)bytes_per_slice_src << endl;
		}

		// DEBUG
		cout << i << '\t' << float(zlib_compressed_summary)/uncompressed_summary << '\t' 
			<< float(rle_compressed_summary)/uncompressed_summary 
			<< "; [" << min_rle_len << ", " << max_rle_len << "]" << endl;

		//throw __FILE__ ": DEBUG";

		fout.write( (char*)dst_buffer, dst_size );

		if(!fout){
			throw __FILE__ ":compress_database_file: Error writing compressed bitslice";
		}
	}

	cerr << "Uncompressed bytes = " << uncompressed_summary << endl;
	cerr << "zlib Compressed bytes = " << zlib_compressed_summary << endl;
	cerr << "RLE Compressed bytes = " << rle_compressed_summary << endl;

	cerr << "zlib Compression ratio = " << double(zlib_compressed_summary)/uncompressed_summary << endl;
	cerr << "RLE Compression ratio = " << double(rle_compressed_summary)/uncompressed_summary << endl;

	// DEBUG
	//cerr << "About to free buffers" << endl;

	// Clear the bitslice buffers
	if(src_buffer != NULL){

		delete [] src_buffer;
		src_buffer = NULL;
	}

	if(dst_buffer != NULL){

		delete [] dst_buffer;
		dst_buffer = NULL;
	}

	if(crc32_src != src_header.crc32){
		throw __FILE__ ":compress_database_file: Invalid CRC32 value for source database";
	}

	// Record the start of the information section
	dst_header.info_start = fout.tellp();

	// We now need to write the filter metadata locations (i.e. the positions of each FilterInfo structure in the
	// file). Since we don't know the sizes yet, write a dummy array as a placeholder that we will
	// fill in later with the actual values.
	unsigned long int *info_loc_buffer = new unsigned long int [dst_header.num_filter];
	
	if(info_loc_buffer == NULL){
		throw __FILE__ ":compress_database_file: Unable to allocate metadata location buffer for destination database file";
	}
	
	fout.write( (char*)info_loc_buffer, dst_header.num_filter*sizeof(unsigned long int) );
		
	if(!fout){
		throw __FILE__ ":compress_database_file: Error writing dummy metadata location buffer to destination databaes file";
	}

	// Skip the info_loc_buffer sections in the input file. We don't need this, as we are reading 
	// the information sequentially.
	fin.seekg( size_t( fin.tellg() ) + sizeof(unsigned long int)*src_header.num_filter);

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// Write the Bloom filter info from the first source file
	for(size_t i = 0;i < src_header.num_filter;++i){

		FilterInfo info;

		binary_read(fin, info);

		// Record the starting file location of this filters metadata in the output file
		info_loc_buffer[i] = fout.tellp();

		binary_write(fout, info);

		if( !fout ){
			throw __FILE__ ":compress_database_file: Error writing Bloom filter info to destination file";
		}
	}
	
	// Rewind to write the compressed sizes and info_loc information for both destination files
	fout.seekp(0);

	binary_write(fout, dst_header);

	if(!fout){
		throw __FILE__ ":compress_database_file: Error writing compressed database file header (final)";
	}

	fout.write( (char*)compressed_bytes, num_bitslice );

	if(!fout){
		throw __FILE__ ":compress_database_file: Error writing compressed database file compressed slice sizes";
	}

	fout.seekp(dst_header.info_start);
	
	fout.write( (char*)info_loc_buffer, dst_header.num_filter*sizeof(unsigned long int) );
		
	if(!fout){
		throw __FILE__ ":compress_database_file: Error writing metadata location buffer to compressed file";
	}

	fout.close();

	if(info_loc_buffer != NULL){

		delete [] info_loc_buffer;
		info_loc_buffer = NULL;
	}

	if(compressed_bytes != NULL){

		delete [] compressed_bytes;
		compressed_bytes = NULL;
	}
}

float bloom_distance(const BitVector &m_a, const BitVector &m_b)
{
	const size_t len = m_a.get_num_bits();

	if(	len != m_b.get_num_bits() ){
		throw __FILE__ ":bloom_distance: Unequal BitVector lengths";
	}

	float ret = 0.0;

	for(size_t i = 0;i < len;++i){
		ret += ( m_a.get_bit(i) == m_b.get_bit(i) );
	}

	// Return the Jaccard distance
	return 1.0 - ret/len;
}

// Return the number of bits in the compressed output
// A return value of 0 indicate no compression (and that not output buffer has been allocated).
// A return value > 0 indicates compression (and an output buffer has been allocated and must be deallocated by the caller).
unsigned int compress_rle(unsigned char** m_out, const unsigned char* m_in, const unsigned int &m_in_bits,
	const unsigned int &m_num_payload_bits, const unsigned int &m_num_count_bits, const unsigned int &m_min_count)
{

	//cerr << "#####################" << endl;

	// Do we have enough bits to compress?
	if(m_in_bits < m_num_payload_bits){
		return 0;
	}

	const unsigned int max_count = (m_min_count - 1) + (1 << m_num_count_bits);

	//cerr << "m_in_bits = " << m_in_bits << endl;
	//cerr << "m_num_count_bits = " << m_num_count_bits << endl;
	//cerr << "m_num_payload_bits = " << m_num_payload_bits << endl;
	//cerr << "m_min_count = " << m_min_count << endl;
	//cerr << "max_count = " << max_count << endl;

	PayloadBuffer mask = 0;

	// Build the mask on the fly
	for(unsigned int i = 0;i < m_num_payload_bits;++i){
		mask |= PayloadBuffer(1) << i;
	}

	if(m_out == NULL){
		throw __FILE__ "compress_rle: Invalid output buffer";
	}

	const unsigned int max_num_output_bytes = m_in_bits/8 + (m_in_bits%8 > 0);

	//cerr << "max_num_output_bytes = " << max_num_output_bytes << endl;

	// The calling function is responsible for freeing this memory
	*m_out = new unsigned char [max_num_output_bytes]; // Max buffer size in bytes

	if(*m_out == NULL){
		throw __FILE__ ":compress_rle: Unable to allocate compressed output buffer";
	}

	// Unset all bits
	memset(*m_out, 0, max_num_output_bytes);

	PayloadBuffer curr_bits = 0;
	PayloadBuffer prev_bits = 0;

	unsigned int prev_num_payload_bits = 0;
	unsigned int curr_num_payload_bits = 0;

	unsigned int out_bit = 0;

	#define	INCREMENT_OUTPUT() \
		if(++out_bit >= m_in_bits){\
			delete [] *m_out;\
			*m_out = NULL;\
			return 0;\
		}

	//cerr << "m_num_payload_bits = " << m_num_payload_bits << endl;

	unsigned int i = 0;

	while(i < m_num_payload_bits){

		prev_bits = (prev_bits << 1) | ( ( m_in[i/8] >> (i%8) ) & 1 );

		++prev_num_payload_bits;
		++i;
	}

	CountBuffer count = 1; // The count is one for the value stored in prev_bits

	// Count the number of input bits that have been successfully stored
	unsigned int num_stored_bits = 0;

	// Format:
	// [is RLE? bit] ? [payload bits][count bits] : [payload bits]
	while(num_stored_bits < m_in_bits){

		//cerr << "raw in bit [" << i << "] = " << ( ( m_in[i/8] >> (i%8) ) & 1 ) << endl;

		//cerr << "num_stored_bits == " << num_stored_bits << endl;
		//cerr << "i = " << i << " out of " << m_in_bits << endl;

		if(i < m_in_bits){

			curr_bits = (curr_bits << 1) | ( ( m_in[i/8] >> (i%8) ) & 1 );
			
			++curr_num_payload_bits;
			++i;
		}

		// Compress in chunks of num_payload_bits (or the remainder if
		// we read the end of the input buffer)
		const bool is_full_payload = (curr_num_payload_bits == m_num_payload_bits);

		//cerr << curr_num_payload_bits << '\t' << m_num_payload_bits << endl;

		if( is_full_payload || (i == m_in_bits) ){

			curr_bits &= mask;

			count += (prev_bits == curr_bits) && is_full_payload;

			//cerr << "i = " << i << ", count = " << count << "; "
			//	<< (is_full_payload ? "full payload" : "partial payload") << endl;
			//cerr << "prev_bits = " << prev_bits << endl;
			//cerr << "curr_bits = " << curr_bits << endl;

			// When (count > max_count), then count == (max_count + 1)
			if( (prev_bits != curr_bits) || (count > max_count) || ( i == m_in_bits ) ){

				if(count >= m_min_count){

					//cerr << "compress: RLE[" << i << "] " << prev_bits << " x " << count << endl;

					//cout << "out_bit (a) = " << out_bit << endl;

					// Set the RLE bit. 
					// Be careful dereferencing, as '[]' binds tighter than '*'.
					(*m_out)[out_bit/8] |= (1 << out_bit%8);
					INCREMENT_OUTPUT();

					//cout << "out_bit (b) = " << out_bit << endl;

					// Set the payload bits one bit at a time (in case we need to cross
					// a byte boundary)
					for(unsigned int j = 0;j < m_num_payload_bits;++j){
						
						// Write the payload bits, least significant bit first.
						// Be careful dereferencing, as '[]' binds tighter than '*'.
						(*m_out)[out_bit/8] |= ( ( (prev_bits >> j) & 1) << out_bit%8);
						INCREMENT_OUTPUT();
					}

					//cout << "out_bit (c) = " << out_bit << endl;

					//cout << "\tcompress RLE count = " << count << endl;

					// Update the stopping criteria
					num_stored_bits += count*m_num_payload_bits;

					//cerr << "\tcompress RLE m_in_bits = " << num_stored_bits << " of " 
					//		<< m_in_bits << endl;

					// Set the count bits one bit at a time (in case we need to cross
					// a byte boundary)
					count -= m_min_count;

					for(unsigned int j = 0;j < m_num_count_bits;++j){
						
						// Write the count bits, least significant bit first.
						// Be careful dereferencing, as '[]' binds tighter than '*'.
						(*m_out)[out_bit/8] |= ( (count >> j) & 1 ) << out_bit%8;
						INCREMENT_OUTPUT();
					}
				}
				else{

					//cerr << "compress: NO RLE (out_bit = " << out_bit << ") " << endl;
					//cerr << "\tcompress no RLE count = " << count << endl;

					//cerr << "NO RLE[" << i << "] " << prev_bits << " x " << count << endl;

					for(unsigned int j = 0;j < count;++j){

						// Do not set the RLE bit (i.e. leave it unset) but increment the
						// outbit to leave the RLE bit == 0
						//*out |= (0 << out_bit);
						INCREMENT_OUTPUT();

						//cout << "out_bit (e) = " << out_bit << endl;

						//cerr << "curr_bits = " << curr_bits << endl;

						for(unsigned int k = 0;k < prev_num_payload_bits;++k){

							//cout << "\tout bit = " << out_bit << "; out_byte = " << out_bit/8 
							//	<< "; m_in_bits = " << m_in_bits << endl;

							// Write the payload bits, least significant bit first.
							// Be careful dereferencing, as '[]' binds tighter than '*'.
							(*m_out)[out_bit/8] |= ( (prev_bits >> k) & 1 ) << out_bit%8;
							INCREMENT_OUTPUT();
						}

						//cout << "out_bit (f) = " << out_bit << endl;

						// Update the stopping criteria
						num_stored_bits += prev_num_payload_bits;

						//cerr << "\tcompress NO RLE m_in_bits = " << num_stored_bits << " of " 
						//	<< m_in_bits << endl;

						if(num_stored_bits == m_in_bits){
							break;
						}
					}
				}

				prev_bits = curr_bits;
				prev_num_payload_bits = curr_num_payload_bits;
				count = 1;

				//cerr << "num_stored_bits = " << num_stored_bits << " out of " << m_in_bits << endl;
			}

			curr_num_payload_bits = 0;
		}
	}

	return out_bit;
}

unsigned int compress_rle(const unsigned char* m_in, const unsigned int &m_len,
	const unsigned int &m_num_payload_bits, const unsigned int &m_num_count_bits, const unsigned int &m_min_count)
{
	unsigned int num_bits = 0;

	const unsigned int max_count = (m_min_count - 1) + (1 << m_num_count_bits);

	PayloadBuffer mask = 0;

	// Build the mask on the fly
	for(unsigned int i = 0;i < m_num_payload_bits;++i){
		mask = (mask << 1) | 1;
	}

	unsigned int prev_bits = 0;
	unsigned int count = 0;

	// Loop over bytes
	unsigned int curr_bits = 0;

	unsigned int bit_index = 0;

	for(unsigned int i = 0;i < m_len;++i){

		// Loop over bits
		for(unsigned int j = 0;j < 8;++j){

			++bit_index;

			curr_bits = (curr_bits << 1) | ( ( m_in[i] >> j ) & 1 );
			
			// Compress in chunks of num_payload_bits
			if(bit_index%m_num_payload_bits == 0){

				curr_bits &= mask;

				if(count == 0){

					prev_bits = curr_bits;
					count = 1;
				}
				else{

					if(prev_bits == curr_bits){
						++count;
					}
					else{
						
						if(count >= m_min_count){
							num_bits += 1 + m_num_payload_bits + m_num_count_bits;
						}
						else{
							num_bits += count*(1 + m_num_payload_bits);
						}

						prev_bits = curr_bits;
						count = 1;
					}

					if(count == max_count){

						num_bits += 1 + m_num_payload_bits + m_num_count_bits;
						count = 0;
					}
				}
			}
		}
	}

	// Include any remainder bits
	if(count > 0){

		if(count >= m_min_count){
			num_bits += 1 + m_num_payload_bits + m_num_count_bits;
		}
		else{
			num_bits += count*(1 + m_num_payload_bits);
		}
	}

	return num_bits/8 + ( (num_bits%8 != 0) ? 1 : 0);
}

void decode_rle(unsigned char* m_out, const unsigned int &m_out_bits, 
	const unsigned char* m_in, 
	const unsigned int &m_num_payload_bits, 
	const unsigned int &m_num_count_bits, 
	const unsigned int &m_min_count)
{
	if( (m_out == NULL) || (m_in == NULL) ){
		throw __FILE__ "decode_rle: Invalid input or output buffer";
	}

	// Set all of the bits in the output buffer to zero
	memset( m_out, 0, m_out_bits/8 + (m_out_bits%8 > 0) );

	unsigned int in_bit = 0;
	unsigned int out_bit = 0;

	//cerr << "-------------------" << endl;

	// Format:
	// [is RLE? bit] ? [payload bits][count bits] : [payload bits]
	while(out_bit < m_out_bits){

		const bool is_rle = ( m_in[in_bit/8] >> in_bit%8 ) & 1;

		++in_bit;

		// Unpack the payload bits, one bit at a time, and store the results
		// in a temporary buffer. Need to add checks to make sure that
		// the maximum number of bits in PayloadBuffer is <= m_num_payload_bits.
		PayloadBuffer buffer = 0;

		for(unsigned int i = 0;i < m_num_payload_bits;++i){

			buffer |= ( ( m_in[in_bit/8] >> in_bit%8 ) & 1 ) << (m_num_payload_bits - 1 - i);

			++in_bit;
		}

		//cerr << "\tpayload = " << buffer << endl;

		if(is_rle){ // RLE encoded

			//cerr << "decode is RLE: " << in_bit << endl;

			// Unpack the count bits, one bit at a time, and store the results
			// in a temporary buffer. Need to add checks to make sure that
			// the maximum number of bits in CountBuffer is <= m_num_count_bits.
			CountBuffer count = 0;

			for(unsigned int i = 0;i < m_num_count_bits;++i){

				count |= ( (m_in[in_bit/8] >> in_bit%8) & 1 ) << i;
				++in_bit;
			}

			// Add the count offset
			count += m_min_count;

			//cout << "\tdecode count = " << count << endl;

			// Write the contents of the payload buffer "count" times to the output buffer
			for(CountBuffer i = 0;i < count;++i){

				for(unsigned int j = 0;j < m_num_payload_bits;++j){

					m_out[out_bit/8] |= ( (buffer >> j) & 1 ) << out_bit%8;
					++out_bit;
				}
			}
		}
		else{ // *Not* RLE encoded

			//cerr << "decode: NOT RLE: " << in_bit << endl;

			// Unpack the payload bits, one bit at a time
			for(unsigned int i = 0;i < m_num_payload_bits;++i){

				//cerr << "\tin bit [" << in_bit << "] = " << ( (buffer >> i) & 1 ) << endl;

				m_out[out_bit/8] |= ( (buffer >> i) & 1 ) << out_bit%8;
				++out_bit;
			}
		}
	}
}