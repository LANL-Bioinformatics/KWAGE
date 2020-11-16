#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>
#include "word.h"
#include "hash.h"
#include "bloom.h"
#include "maestro.h"
#include <ncbi-vdb/NGS.hpp> // For openReadCollection

// For non-temporal memory access that bypasses the cache
#include <smmintrin.h>
#include <emmintrin.h>

using namespace std;

int main(int argc, char *argv[])
{
	try{

		time_t profile = time(NULL);

		const size_t num_bloom = 257;

		deque<string> bloom_files;

		for(size_t i = 0;i < num_bloom;++i){
			
			stringstream filename;

			filename << "temp/filter." << i << ".bloom";

			bloom_files.push_back( filename.str() );
		}

		BloomParam param;

		param.hash_func = MURMUR_HASH_32;
		param.log_2_filter_len = 18;
		param.kmer_len = 31;
		param.num_hash = 3;

		// Create Bloom filters
		if(false){

			for(deque<string>::const_iterator i = bloom_files.begin();
				i != bloom_files.end();++i){

				cerr << "Writing filter " << *i << endl;

				ofstream fout(i->c_str(), ios::binary);

				if(!fout){
					throw __FILE__ ":main: Unable to open output filter for writing";
				}

				const size_t num_bits = 1ULL << param.log_2_filter_len;

				BloomFilter filter(param);

				filter.unset_all_bits();

				for(size_t j = 0;j < num_bits;++j){

					if(rand() % 2 == 1){
						filter.set_bit(j);
					}
				}

				// Updat the crc32 value
				filter.update_crc32();

				// Don't bother with any metadata
				binary_write(fout, filter);
			}
		}

		const bool ret = build_db("temp/debug_new.db", param, bloom_files);

		if(ret){
			cerr << "Successfully wrote the database file" << endl;
		}
		else{
			cerr << "*Failed* to write the database file" << endl;
		}

		profile = time(NULL) - profile;

		cerr << "Complete in " << profile << " sec" << endl;
	}
	catch(const char* error){

		cerr << "Caught the error: " << error << endl;
		return EXIT_FAILURE;
	}
	catch(...){

		cerr << "Caught an unhandled error" << endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
