// Check the integrity of an input Bloom filter and optionally 
// delete the source Bloom filter data.
// J. D. Gans
// Bioscience Division, B-11
// Los Alamos National Laboratory
// Tue Jan 28 15:41:53 2020

#include <iostream>
#include <algorithm>
#include <iomanip>

#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

#include "file_util.h"
#include "bloom.h"
#include "bigsi++.h"

using namespace std;

int main(int argc, char *argv[])
{
	try{
		
		// Command line arguments
		// --bloom <input Bloom filter name>
		// [--clean <SRA data directory to delete>]
		// [-v (verbose output)]

		const char* options = "vh?";
		int config_opt = 0;
		int long_index = 0;

		struct option long_opts[] = {
			{"bloom", true, &config_opt, 1},
			{"clean", true, &config_opt, 2},
			{0,0,0,0} // Terminate options list
		};

		int opt_code;
		opterr = 0;
		
		bool print_usage = (argc == 1);
		
		string bloom_filename;
		string sra_data_dir;
		bool verbose = false;
		
		while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

			switch( opt_code ){
				case 0:
					if(config_opt == 1){ // --bloom 

						bloom_filename = optarg;
						break;
					}

					if(config_opt == 2){ // --clean 

						sra_data_dir = optarg;
						break;
					}

					cerr << "Unknown flag!" << endl;
					break;
				case 'v':
					verbose = true;
					break;
				case '?':
				case 'h':
					print_usage = true;
					break;
				default:
					cerr << '\"' << (char)opt_code << "\" is not a valid option!" << endl;
					break;
			};
		}

		if(print_usage){

			cerr << "Usage for checkbloom (v. " << CHECK_BLOOM_VERSION << "):" << endl;
			cerr << "\t--bloom <input Bloom filter name>" << endl;
			cerr << "\t[--clean <SRA data directory to delete>]" << endl;
			cerr << "\t[-v (verbose output)]" << endl;

			return EXIT_SUCCESS;
		}
		
		if( bloom_filename.empty() ){
			
			cerr << "Please specify an input Bloom filter file (--bloom)" << endl;
			return EXIT_SUCCESS;
		}

		if( !sra_data_dir.empty() && !is_dir(sra_data_dir) ){
		
			cerr << "Please specify a valid directory to clean (--clean)" << endl;
			return EXIT_FAILURE;
		}
		
		ifstream fin(bloom_filename.c_str(), ios::binary);
		
		if(!fin){
			
			cerr << "Unable to open Bloom filter file for reading: " << bloom_filename << endl;
			return EXIT_FAILURE;
		}
		
		BloomFilter filter;
		
		binary_read(fin, filter);
				
		if( filter.empty() && filter.get_param().empty() ){
			
			if(verbose){
				cerr << "Empty Bloom filter file (unable to satisfy false positive rate)" << endl;
			}
		}
		else{
			if( !filter.test_crc32() ){

				if(verbose){
					cerr << "Checksum failure for: " << bloom_filename << endl;
				}

				// Exit before removing any SRA data!
				return EXIT_FAILURE;
			}
			else{
				if(verbose){
					cerr << "Checksum is valid : " 
						<< std::hex << filter.get_crc32() 
						<< std::dec << endl;
					
					const size_t filter_count = filter.count();
					
					cerr << "Number of set Bloom filter bits = " << filter_count
						<< " (" << (100.0*filter_count)/filter.get_param().filter_len() 
						<< "%)" << endl;
				}
			}
		}
		
		// If the user does not specify a directory to clean up, we can return now
		if( sra_data_dir.empty() ){
			return EXIT_SUCCESS;
		}
		
		FindFiles ff(sra_data_dir);
		
		deque<string> del_dir;

		while( ff.next() ){

			if(unlink( ff.name().c_str() ) != 0){
				
				if(verbose){
					cerr << "Warning: Unable to unlink file: " << ff.name() << endl;
				}
			}
			else{
				if(verbose){
					cerr << "Removed file: " << ff.name() << endl;
				}
			}

			del_dir.push_back( ff.dir() );
		}

		del_dir.push_back(sra_data_dir);

		// Make the list of directories to delete unique
		sort( del_dir.begin(), del_dir.end() );
		del_dir.erase( unique( del_dir.begin(), del_dir.end() ), del_dir.end() );

		sort( del_dir.begin(), del_dir.end(), sort_by_length() );

		for(deque<string>::const_reverse_iterator i = del_dir.rbegin();i != del_dir.rend();++i){

			if(rmdir( i->c_str() ) != 0){

				if(verbose){
					cerr << "Warning: Error removing sra run directory -> " 
						<< strerror(errno) << endl;
				}
			}
			else{
				if(verbose){
					cerr << "Removed directory: " << *i << endl;
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
