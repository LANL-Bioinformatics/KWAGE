#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <zlib.h>

#include "bloom.h"
#include "binary_io.h"

using namespace std;

inline unsigned char pop_count(unsigned char m_byte)
{
	return 
		(m_byte & 1) + 
		( (m_byte >> 1) & 1) +
		( (m_byte >> 2) & 1) +
		( (m_byte >> 3) & 1) +
		( (m_byte >> 4) & 1) +
		( (m_byte >> 5) & 1) +
		( (m_byte >> 6) & 1) +
		( (m_byte >> 7) & 1);
}

int main(int argc, char *argv[])
{
	try{
		if(argc != 3){

			cerr << "Usage: " << argv[0] << " <Bloom filter file 1> <Bloom filter file 2>" << endl;
			return EXIT_SUCCESS;
		}

		ifstream fin_1(argv[1], ios::binary);

		if(!fin_1){

			cerr << "Unable to open file 1: " << argv[1] << endl;
			return EXIT_FAILURE;
		}

		ifstream fin_2(argv[2], ios::binary);

		if(!fin_2){

			cerr << "Unable to open file 2: " << argv[2] << endl;
			return EXIT_FAILURE;
		}

		// Make sure that all of the Bloom filters have been completely written
		unsigned char status;

		binary_read(fin_1, status);

		if(status != BLOOM_MAGIC_COMPLETE){

			cerr << "Bloom filter 1 (" << argv[1] << ") is not complete!" << endl;
			return EXIT_FAILURE;
		}

		binary_read(fin_2, status);

		if(status != BLOOM_MAGIC_COMPLETE){

			cerr << "Bloom filter 2 (" << argv[2] << ") is not complete!" << endl;
			return EXIT_FAILURE;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read the BloomParam data from each filter. We won't be storing the this information for each filter,
		// but do need to make sure that the filter paramters are all the *same*.

		BloomParam param1;
		BloomParam param2;

		binary_read(fin_1, param1);

		if(!fin_1){

			cerr << "Error reading parameters from Bloom filter 1 (" << argv[1] << ")" << endl;
			return EXIT_FAILURE;
		}

		binary_read(fin_2, param2);

		if(!fin_2){

			cerr << "Error reading parameters from Bloom filter 2 (" << argv[2] << ")" << endl;
			return EXIT_FAILURE;
		}

		if(param1 != param2){

			cerr << "Inconsistent Bloom filter parameters" << endl;
			cerr << "kmer_len = " << param1.kmer_len << " for 1;  " << param2.kmer_len << " for 2" << endl;
			cerr << "log_2_filter_len = " << param1.log_2_filter_len << " for 1; expected " << param2.log_2_filter_len << " for 2" << endl;
			cerr << "num_hash = " << param1.num_hash << " for 1; expected " << param2.num_hash << " for 2" << endl;
			cerr << "hash_func = " << hash_name(param1.hash_func) << " for 1; expected " << hash_name(param2.hash_func) << " for 2" << endl;

			return EXIT_FAILURE;
		}

		const size_t filter_len = param1.filter_len();

		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read the CRC data from each filter. These crc32 values will be used to validate the Blooom filter
		// data that we read. After all Bloom filters have been read, the reference and computed CRC32 values
		// should be equal. Note that before we compute the running crc32 values, we must initialize
		// the crc32 value with a call to crc32_z(0L, Z_NULL, 0).
		unsigned int crc32_file_1;
		unsigned int crc32_file_2;

		binary_read(fin_1, crc32_file_1);

		if(!fin_1){

			cerr << "Unable to read crc32 value from Bloom filter 1 (" << argv[1] << ")" << endl;
			return EXIT_FAILURE;
		}

		binary_read(fin_2, crc32_file_2);

		if(!fin_2){

			cerr << "Unable to read crc32 value from Bloom filter 2 (" << argv[2] << ")" << endl;
			return EXIT_FAILURE;
		}

		if(crc32_file_1 == crc32_file_2){

			cerr << "The crc32 values are the same for both Bloom filters (" 
				<< std::hex << crc32_file_1 << std::dec << ")" << endl;
		}
		else{

			cerr << "The Bloom filters have different crc32 values" << endl;
			cerr << "\tBloom filter 1 (" << argv[1] << ") crc32 = " << std::hex << crc32_file_1 << std::dec << endl;
			cerr << "\tBloom filter 2 (" << argv[2] << ") crc32 = " << std::hex << crc32_file_2 << std::dec << endl;
		}

		unsigned int crc32_computed_1 = crc32_z(0L, Z_NULL, 0);
		unsigned int crc32_computed_2 = crc32_z(0L, Z_NULL, 0);

		FilterInfo info1;
		FilterInfo info2;

		binary_read(fin_1, info1);

		if(!fin_1){

			cerr << "Unable to read the FilterInfo from Bloom filter 1 (" << argv[1] << ")" << endl;
			return EXIT_FAILURE;
		}

		binary_read(fin_2, info2);

		if(!fin_2){

			cerr << "Unable to read the FilterInfo from Bloom filter 2 (" << argv[2] << ")" << endl;
			return EXIT_FAILURE;
		}

		const size_t buffer_size = 1048576; // 256 MB buffer

		if(buffer_size%BitVector::BITS_PER_BLOCK != 0){
			throw __FILE__ ":main: Buffer size must be a multiple of BitVector::BITS_PER_BLOCK";
		}

		unsigned char *buffer1 = new unsigned char [buffer_size];

		if(buffer1 == NULL){
			throw __FILE__ ":main: Unable to allocate buffer 1";
		}

		unsigned char *buffer2 = new unsigned char [buffer_size];

		if(buffer2 == NULL){
			throw __FILE__ ":main: Unable to allocate buffer 2";
		}

		size_t diff_bits = 0;

		for(size_t i = 0;i < filter_len;i += buffer_size){
			
			const size_t num_byte = min(buffer_size, filter_len - i)/BitVector::BITS_PER_BLOCK;

			fin_1.read( (char*)buffer1, num_byte );

			if(!fin_1){
				throw __FILE__ ":main: Error reading filter 1 bytes";
			}

			// Keep a running update of the per-filter CRC32 values
			crc32_computed_1 = ::crc32_z(crc32_computed_1, buffer1, num_byte);

			fin_2.read( (char*)buffer2, num_byte );

			if(!fin_2){
				throw __FILE__ ":main: Error reading filter 2 bytes";
			}

			// Keep a running update of the per-filter CRC32 values
			crc32_computed_2 = ::crc32_z(crc32_computed_2, buffer2, num_byte );

			for(size_t j = 0;j < num_byte;++j){

				// Count the bits that *disagree*
				diff_bits += pop_count(buffer1[j] ^ buffer2[j]);
			}
		}

		cerr << "The Bloom filters differ by " << diff_bits << " bits of out " << filter_len << " bits: " 
			<< (100.0*diff_bits)/filter_len << "%" << endl;

		if(crc32_computed_1 != crc32_file_1){

			cerr << "The crc32 disagreement for Bloom filter 1: " << argv[1] << endl;
			cerr << "\tComputed crc32: " << std::hex << crc32_computed_1 << std::dec << endl;
			cerr << "\tFile crc32: " << std::hex << crc32_file_1 << std::dec << endl;
		}

		if(crc32_computed_2 != crc32_file_2){

			cerr << "The crc32 disagreement for Bloom filter 2: " << argv[1] << endl;
			cerr << "\tComputed crc32: " << std::hex << crc32_computed_2 << std::dec << endl;
			cerr << "\tFile crc32: " << std::hex << crc32_file_2 << std::dec << endl;
		}

		if(buffer1 != NULL){

			delete [] buffer1;
			buffer1 = NULL;
		}

		if(buffer2 != NULL){

			delete [] buffer2;
			buffer2 = NULL;
		}
	}
	catch(const char *error){

		cerr << "Caught the error: " << error << endl;
		return EXIT_FAILURE;
	}
	catch(...){

		cerr << "Caught an unhandled error!" << endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
