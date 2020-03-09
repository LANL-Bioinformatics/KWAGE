// Build a BIGSI Bloom filter database from input Bloom filter files
// J. D. Gans
// Bioscience Division, B-10
// Los Alamos National Laboratory
// Fri Oct 25 10:22:49 2019

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdlib.h>
#include <signal.h>
#include <math.h>

#include "bigsi++.h"
#include "options.h"
#include "parse_sequence.h"
#include "bloom.h"
#include "hash.h"
#include "sort.h"
#include "mpi_util.h"
#include "file_util.h"
#include "update.h"
#include "database.h"
#include "keys.h"

using namespace std;

// Global variables for MPI
int mpi_numtasks;
int mpi_rank;

double start_time;

void terminate_program(int m_sig);
string report_run_time();

int main(int argc, char *argv[])
{
	try{

		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_numtasks);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
                
		signal( SIGINT, terminate_program );
        	signal( SIGTERM, terminate_program );
        	signal( SIGSEGV, terminate_program );
                
        	start_time = MPI_Wtime();
		
		BuildOptions opt;
		
		// Only rank 0 reads the command line arguments
		if(mpi_rank == 0){
			opt.load(argc, argv);
		}
		
		broadcast(opt, mpi_rank, 0);
		
		if(opt.quit){
			
			MPI_Finalize();
			return EXIT_SUCCESS;
		}
		
		// Only rank 0 logs progress information
		UpdateInfo progress(mpi_rank == 0);
			
		if(mpi_rank == 0){
			
			if( !opt.log_filename.empty() ){
			
				try{
					progress.log(opt.log_filename);
				}
				catch(...)
				{
					cerr << "Unable to open logfile: " << opt.log_filename << endl;
					opt.quit = true;
				}
			}
		
			// Prepare the output directory structure by making sure that the
			// root-level directory exists:
			//	root/
			//		k_hashname/
			//			ceil( log_2(M1) )/
			//			ceil( log_2(M2) )/
			//			ceil( log_2(M3) )/
			//				num_hash_func
			//			...
			if( !make_dir(opt.output_dir) ){
				
				cerr << "Unable to create the output directory: " << opt.output_dir << endl;
				opt.quit = true;
			}
			
			#ifdef NOT_NOW
			stringstream ssout;
			
			ssout << opt.output_dir << PATH_SEPARATOR 
				<< opt.kmer_len << "_" << hash_name(opt.hash_func);
			
			opt.output_dir = ssout.str();
			
			if( !make_dir(opt.output_dir) ){
				
				cerr << "Unable to create the k-mer & hash output sub-directory: " << opt.output_dir << endl;
				opt.quit = true;
			}
			#endif // NOT_NOW
			
		}
		
		broadcast(opt.quit, mpi_rank, 0);
		
		if(opt.quit){
			
			MPI_Finalize();
			return EXIT_SUCCESS;
		}
		
		// Let all ranks know the new output sub-directory
		broadcast(opt.output_dir, mpi_rank, 0);
		
		const size_t num_filter = opt.filter_input.size();
		
		// Find the Bloom filter parameters for each input Bloom filter.
		unordered_map< BloomParam, deque<size_t> > param_to_file_index;
		size_t num_empty_filter = 0;
		
		if(mpi_rank == 0){
		
			const size_t update_every = max( size_t(1), num_filter/100 );
			
			for(size_t i = 0;i < num_filter;++i){

				// Read the Bloom filter parameters associated with
				// the selected file
				ifstream fin(opt.filter_input[i].c_str(), ios::binary);

				if(!fin){
					throw __FILE__ ":main: Unable to open Bloom filter file for parameter extraction";
				}

				BloomParam param;

				binary_read(fin, param);

				// Some of the Bloom filter files may be empty if the source SRA files
				// contained too many kmers for the allowed maximum Bloom filter length.
				// Don't include these empty files in the final database (but do report
				// the number of files that were skipped)
				if( param.empty() ){
					++num_empty_filter;
				}
				else{
					param_to_file_index[param].push_back(i);
				}
				
				if(i%update_every == 0){
					
					progress << "Scanning " << opt.filter_input[i] 
						<< " (" << (100.0*i)/num_filter << "%)" ;
					progress.flush();
				}
			}
			
			progress << "Found " << param_to_file_index.size() 
				<< " unique parameter sets for the " << num_filter << " input Bloom filters";
			progress.flush();
			
			progress << "Skipped " << num_empty_filter
				<< " empty Bloom filters";
			progress.flush();
			
			progress.close();
		}
		
		broadcast(param_to_file_index, mpi_rank, 0);
		
		// Make sure that all ranks agree on the parameter order
		const vector<BloomParam> param = keys(param_to_file_index);
		
		for(vector<BloomParam>::const_iterator p = param.begin();p != param.end();++p){
			
			progress << "Building database for the \"" << hash_name(p->hash_func) 
				<< "\" hash; log2(len) = " << p->log_2_filter_len << "; "
				<< p->num_hash << " hash function" << ( (p->num_hash > 1) ? "s" : "");
			progress.flush();
			progress.close();
			
			unordered_map< BloomParam, deque<size_t> >::const_iterator iter = param_to_file_index.find(*p);
			
			if( iter == param_to_file_index.end() ){
				throw __FILE__ ":main: Unable to lookup Bloom filter parameters";
			}
			
			const size_t num_local_filter = iter->second.size();
			
			// Each rank loads a unique subset of Bloom filters. Recall that we are using MPI primarily
			// to distribute the memory across multiple ranks
			deque<BloomFilter> db;
			
			for(size_t i = 0;i < num_local_filter;++i){
				
				if( i%( size_t(mpi_numtasks) ) == size_t(mpi_rank) ){
					
					db.push_back( BloomFilter() );
					
					ifstream fin(opt.filter_input[ iter->second[i] ].c_str(), ios::binary);
					
					binary_read( fin, db.back() );
				}
			}
			
			// Merge these Bloom filters with any existing database records
			merge_bloom_filters(db, *p, progress, opt);
		}
		
		if(mpi_rank == 0){
			cerr << report_run_time() << endl;
		}
		
		MPI_Finalize();
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

// Report our rank, the signal we caught and the time spent running
void terminate_program(int m_sig)
{
        cerr << "[" << mpi_rank << "] caught signal " << m_sig << endl;
        cerr << report_run_time() << endl;
        
        MPI_Abort(MPI_COMM_WORLD, 0);
}

// Run time computes the total run time. The results are formatted as a string.
string report_run_time()
{
        double elapsed_time = MPI_Wtime() - start_time; // In sec
        
        const double elapsed_sec = fmod(elapsed_time, 60.0);
        
        elapsed_time = (elapsed_time - elapsed_sec)/60.0; // In min
        
        const double elapsed_min = fmod(elapsed_time, 60.0);
        elapsed_time = (elapsed_time - elapsed_min)/60.0; // In hour
        
        const double elapsed_hour = fmod(elapsed_time, 24.0);
        elapsed_time = (elapsed_time - elapsed_hour)/24.0; // In day
        
        stringstream sout;
        
        sout << "Run time is " 
                << elapsed_time 
                << " days, "
                << elapsed_hour
                << " hours, "
                << elapsed_min
                << " min and "
                << elapsed_sec
                << " sec";
        
        return sout.str();
}
