// Convert SRA sequence files and metadata into Bloom filters
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <deque>
#include <unordered_map>
#include <algorithm>

#include <stdlib.h>
#include <getopt.h>
#include <signal.h>
#include <math.h>

#include <omp.h>
#include <ncbi-vdb/NGS.hpp> // For openReadCollection

#include "options.h"
#include "bloom.h"
#include "hash.h"
#include "word.h"
#include "parse_sequence.h"
#include "file_util.h"
#include "sort.h"

using namespace std;

// Global variables for MPI
int mpi_numtasks;
int mpi_rank;

double start_time;

void terminate_program(int m_sig);
string report_run_time();

// This hash function is only used to uniquely assign kmers to different
// mpi ranks. Is it *not* used to compute Bloom filters.
inline int kmer_hash(Word m_kmer, const int &m_num_task)
{
	// From http://burtleburtle.net/bob/hash/integer.html
	// This website has a number of pretty good hash functions.
	// Even though these functions expect 32 bit input, they still seem
	// to work well for the (up to) 64 bit inputs used below (I've currently
	// only tested with 48 bit Words).
	// Hash #1
	//w = (w^0xdeadbeef) + (w<<4);
	//w = w ^ (w>>10);
	//w = w + (w<<7);
	//w = w ^ (w>>13);
	//return w % m_num_tasks;

	// Hash #2
	//w = w ^ (w>>4);
    	//w = (w^0xdeadbeef) + (w<<5);
    	//w = w ^ (w>>11);
    	//return w % m_num_tasks;

	// Hash #3
    	m_kmer -= (m_kmer << 6);
    	m_kmer ^= (m_kmer >> 17);
    	m_kmer -= (m_kmer << 9);
    	m_kmer ^= (m_kmer << 4);
    	m_kmer -= (m_kmer << 3);
    	m_kmer ^= (m_kmer << 10);
    	m_kmer ^= (m_kmer >> 15);
	
    	return m_kmer % m_num_task;
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_numtasks);
    	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
                
    	signal( SIGINT, terminate_program );
    	signal( SIGTERM, terminate_program );
    	signal( SIGSEGV, terminate_program );
                
    	start_time = MPI_Wtime();
		
	try{
		
		BloomerOptions opt;
		
		// Only rank 0 reads the command line arguments
		if(mpi_rank == 0){
			opt.load(argc, argv);
		}
		
		broadcast(opt, mpi_rank, 0);
		
		if(opt.quit){
			
			MPI_Finalize();
			return EXIT_SUCCESS;
		}
		
		time_t profile = time(NULL);
		
		if( opt.verbose && (mpi_rank == 0) ){
			
			cout << "K-mer length = " << opt.kmer_len << endl;
			cout << "False positive probability = " << opt.false_positive_probability << endl;
			cout << "Minimum k-mer count = " << opt.min_kmer_count << endl;
			cout << "Using the \"" << hash_name(opt.hash_func) << "\" hash function" << endl;
			
			if(opt.save_sra){
				cout << "SRA files will *not* be deleted after Bloomer filter construction!" << endl;
			}
		}
		
		// Find all of the SRA files that are associated with the provided run name
		FindFiles ff(opt.input_dir);

		deque<string> sra_files;
		deque<string> meta_files;

		while( ff.next() ){

			if( find_file_extension(ff.name(), ".sra") ) {
				sra_files.push_back( ff.name() );
			}
			
			if( find_file_extension(ff.name(), ".info") ) {
				meta_files.push_back( ff.name() );
			}
		}

		if( sra_files.empty() ){
			
			cerr << "Did not find any SRA files in " << opt.input_dir << endl;
			
			MPI_Finalize();
			return EXIT_SUCCESS;
		}
		
		if( opt.read_meta_data && (meta_files.size() != 1) ){
			
			cerr << "Did not find any metadata files in " << opt.input_dir << endl;
			
			MPI_Finalize();
			return EXIT_SUCCESS;
		}
		
		// Digest the current sequence group. Trial and error testing on beagle.lanl.gov
		// indicates that 5 reading threads provides the best performance. Binding to 
		// sock (with mpirun --bind-to socket) provides better read performance at the
		// expense of worse parallel sorting performance (so don't do it). 
		
		if( opt.verbose && (mpi_rank == 0) ){
			cout << "Using " << opt.num_file_slice << " threads to read SRA data" << endl;
		}
		
		FilterInfo info;
		
		if(opt.read_meta_data && (mpi_rank == 0) ){
		
			// Only rank 0 needs to read the metadata
			ifstream fin(meta_files[0].c_str(), ios::binary);

			if(!fin){
				throw __FILE__ ":main: Unable to open metadata information file for reading";
			}

			binary_read(fin, info);
		}
		
		vector< deque<Word> > curr_kmers(opt.num_file_slice);

		time_t kmer_profile = time(NULL);
		
		for(deque<string>::const_iterator i = sra_files.begin();i != sra_files.end();++i){

			if( opt.verbose && (mpi_rank == 0) ){
				cout << "Reading kmers from: " << *i << endl;
			}
			
			// For SRA records with additional information (like reference sequences and
			// .sra.vdbcache file), we need to be in a directory that is within
			// two leves of the actual sra file.
			// For a simple test case:
			// 	./foo SRR9008316.sra --> success
			// 	./foo SRR9008316/SRR9008316.sra --> success
			// 	./foo SRR9008316_run/SRR9008316/SRR9008316.sra --> *failure*
			//
			// This behaviour seems like a bug to me, but we can work around it by
			// extracting the path from each SRA filename and changing the working
			// directory to this path.
			
			char* orig_dir = getcwd(NULL, 0);

			if(orig_dir == NULL){
				throw __FILE__ ":main: Unable to obtain original working directory";
			}

			const pair<string /*dir*/, string /*file*/> dir_and_file = split_dir_and_file(*i);
			
			chdir( dir_and_file.first.c_str() );
			
			// Reading multiple chunks in parallel appears to be 
			// slightly faster when the ngs::ReadCollection variable
			// is thread local.		
			//ngs::ReadCollection run(  ncbi::NGS::openReadCollection(dir_and_file.second) );
			//const size_t num_read = run.getReadCount(ngs::Read::all);
			
			#pragma omp parallel num_threads(opt.num_file_slice)
			{
				const size_t num_thread = omp_get_num_threads();
				const size_t tid = omp_get_thread_num();
				
				deque<Word>& local_kmers = curr_kmers[tid];
				
				ngs::ReadCollection run(  ncbi::NGS::openReadCollection(dir_and_file.second) );
				
				// Note that num_read is the number of either paired or
				// unpaired reads. For paired reads, this is half the
				// the number of sequences!
				const size_t num_read = run.getReadCount(ngs::Read::all);
				
				size_t chunk = num_read/num_thread;
											
				// Each thread is assigned an non-overlapping slice of the SRA
				// file to read. Since getReadRange is 1's based, we must add 1 to start.
				const size_t start = 1 + chunk*tid;
				
				if( tid == (num_thread - 1) ){

					// Read any remainder reads with the last thread
					chunk += num_read%num_thread;
				}
								
				ngs::ReadIterator run_iter = 
					ngs::ReadIterator( run.getReadRange ( start, chunk, ngs::Read::all ) );
				
				size_t seq_count = 0;
				size_t read_count = 0;
							
				while( run_iter.nextRead() ){
					
					++read_count;
					
					while( run_iter.nextFragment() ){
						
						++seq_count;
						
						// DEBUG
						if( opt.verbose && (tid == 0) && (mpi_rank == 0) && (seq_count%1000000 == 0) ){
							cout << ((double)(seq_count))/(time(NULL) - kmer_profile) 
								<< " read/sec (" << (100.0*read_count)/chunk << "%)" << endl;
						}
						
						const string seq = run_iter.getFragmentBases().toString();
						
						ForEachDuplexWord(seq, opt.kmer_len)

							if(ValidWord){

								const Word w = CanonicalWord;

								// Partition the words so that each MPI rank gets a 
								// unique subset of words. Note that the kmer_hash uses
								// a different hashing function from the hash used to 
								// construct the Bloom filters
								if(kmer_hash(w, mpi_numtasks) == mpi_rank){
									local_kmers.push_back(w);
								}
							}

						EndWord
					}
				}				
			}
				
			// Return to the original directory after reading each SRA file
			chdir(orig_dir);
		}
		
		kmer_profile = time(NULL) - kmer_profile;
		
		if( opt.verbose && (mpi_rank == 0) ){
			
			cout << "\tDone reading in " << kmer_profile << " sec" << endl;
			cout << "\tFinding abundant kmers ... ";
		}
		
		time_t abundant_profile = time(NULL);
		
		deque<Word> valid_kmers; // kmers that have passed the minimum count threshold
			
		// Count the occurance of each k-mer for frequency-based
		// k-mer filtering (to remove sequencing errors)
		find_abundant_kmers(valid_kmers, curr_kmers, opt.min_kmer_count);
		
		abundant_profile = time(NULL) - abundant_profile;
		
		if( opt.verbose && (mpi_rank == 0) ){
			cout << "done (in " << abundant_profile << " sec)" << endl;
		}
		
		// Discard any kmers remaining in curr_kmers, they did not appear
		// at least opt.min_kmer_count times
		curr_kmers.clear(); // Free memory
		
		unsigned long long int subset_valid_unique_kmer = valid_kmers.size();
		
		// Recall that each rank only has a subset of the unique kmers.
		// Sum all of the unique and valid kmers across all ranks to
		// get the actual count.
		unsigned long long int num_valid_unique_kmer = 0;
		
		MPI_Allreduce(&subset_valid_unique_kmer, &num_valid_unique_kmer, 1, MPI_UNSIGNED_LONG_LONG, 
			MPI_SUM, MPI_COMM_WORLD);
		
		BloomParam param;
		
		try{
			param = optimal_bloom_param(opt.kmer_len,
				num_valid_unique_kmer, 
				opt.false_positive_probability,
				opt.hash_func);
		}
		catch(...){

			// We were unable to find Bloom filter parameters that
			// satisfied the requested false_positive_probability.
			// What should we do? 
			//
			// For now, we will write a special, *empty* Bloom filter file
			// that will indicate that we could not satisfy the requested
			// error rate without exceeding the maximum allowed Bloom filter
			// file size.
			if(mpi_rank == 0){
				
				if(opt.verbose){
					cout << "Unable to satisfy the requested false positive rate. Giving up!" << endl;
				}
				
				BloomFilter filter; // An empty filter
				
				filter.set_info(info); // Set the metadata (so we can inventory failed SRA records)
				
				ofstream fout(opt.output_file, ios::binary);

				if(!fout){
					throw __FILE__ ":main: Unable to open Bloom filter file for writing";
				}

				binary_write(fout, filter);
				
				if( !opt.save_sra ){
				
					// Delete the input directory and all of its contents
					// Remove any associated SRA files (including *.sra and *.sra.vdbcache)
					// and metadata files
					remove_all(opt.input_dir);

					if(opt.verbose){
						cout << "Deleted SRA input data in " << opt.input_dir << endl;
					}
				}
			}
			
			MPI_Finalize();
			return EXIT_SUCCESS;
		}

		const uint32_t filter_len = param.filter_len();
		
		if( opt.verbose && (mpi_rank == 0) ){
			
			cout << "Found " << num_valid_unique_kmer << " unique k-mers" << endl;
			cout << "Bloom filter length = 2^" << param.log_2_filter_len << " = " << filter_len  << endl;
			cout << "Number of hash functions = " << param.num_hash << endl;
		}
		
		// Each rank accumulates the bits for the subset of kmers belonging to the rank.
		// These bit vectors will then be reduced on rank 0.
		BloomFilter filter(param);
		
		filter.unset_all_bits();
		
		#pragma omp parallel
		{
			// A thread local Bloom filter
			BloomFilter local(param);

			local.unset_all_bits();

			#pragma omp for
			for(unsigned long long int i = 0;i < subset_valid_unique_kmer;++i){

				// The BIGSI python implementation starts from 0 when seeding the
				// hash function
				for(size_t h = 0;h < param.num_hash;++h){

					// ** Check hash value against python version -- potential off by one
					// issue involving mod of negative number !! ***
					local.set_bit( bigsi_hash(valid_kmers[i], opt.kmer_len, h, filter_len, opt.hash_func) );
				}
			}

			// The Bloom filter bitwise OR operator copy data if the destination is not the same size as
			// the right hand side.
			#pragma omp critical
			filter |= local;
		}
		
		// Free up space as we go!
		valid_kmers.clear();
		
		// Reduce all of the filter bits to rank 0 using the binary OR operation
		if(mpi_rank == 0){
			MPI_Reduce(MPI_IN_PLACE, filter.ptr(), filter.num_block()*sizeof(BitVector::BLOCK), 
				MPI_BYTE, MPI_BOR, 0, MPI_COMM_WORLD);
		}
		else{
			MPI_Reduce(filter.ptr(), NULL, filter.num_block()*sizeof(BitVector::BLOCK), 
				MPI_BYTE, MPI_BOR, 0, MPI_COMM_WORLD);
		}
		
		// Rank 0 is responsible for writing the Bloom filter
		if(mpi_rank == 0){
		
			// Compute the checksum value to safeguard the Bloom filter data only (metadata and Bloom filter
			// parameters are not included in the crc32).
			const unsigned int checksum = filter.update_crc32();
			
			if(opt.verbose){
				cout << "Bloom filter checksum = " << std::hex << checksum << std::dec << endl;
				cout << "Bloom filter occupancy = " << ( (float)filter.count() )/filter_len << endl;
			}
			
			filter.set_info(info); // Set the metadata
			
			if(opt.verbose){
				
				cout << "Meta data:" << endl;
				cout << "\trun_accession: " << info.run_accession << endl;
				cout << "\texperiment_accession: " << info.experiment_accession << endl;
				cout << "\texperiment_title: " << info.experiment_title << endl;
				cout << "\texperiment_design_description: " << info.experiment_design_description << endl;
				cout << "\texperiment_library_name: " << info.experiment_library_name << endl;
				cout << "\texperiment_library_strategy: " << info.experiment_library_strategy << endl;
				cout << "\texperiment_library_source: " << info.experiment_library_source << endl;
				cout << "\texperiment_library_selection: " << info.experiment_library_selection << endl;
				cout << "\texperiment_instrument_model: " << info.experiment_instrument_model << endl;
				cout << "\tsample_accession: " << info.sample_accession << endl;
				cout << "\tsample_taxa: " << info.sample_taxa << endl;
				
				if( !info.sample_attributes.empty() ){
				
					cout << "\tsample attributes:" << endl;
					
					for(MULTIMAP<string, string>::const_iterator j = info.sample_attributes.begin();
						j != info.sample_attributes.end();++j){
						
						cout << "\t\t" << j->first << ": " << j->second << endl;
					}
				}
				
				cout << "\tstudy_accession: " << info.study_accession << endl;
				cout << "\tstudy_title: " << info.study_title << endl;
				cout << "\tstudy_abstract: " << info.study_abstract << endl;
			}
			
			ofstream fout(opt.output_file, ios::binary);
			
			if(!fout){
				throw __FILE__ ":main: Unable to open Bloom filter file for writing";
			}
			
			binary_write(fout, filter);
			
			// Make sure that this file is closed and commited to disk before we attempt to 
			// delete the source SRA data.
			fout.close();

			if( !opt.save_sra ){
				
				// Delete the input directory and all of its contents
				// Remove any associated SRA files (including *.sra and *.sra.vdbcache)
				// and metadata files
				remove_all(opt.input_dir);
				
				if(opt.verbose){
					cout << "Deleted SRA input data in " << opt.input_dir << endl;
				}
			}
		}
		
		profile = time(NULL) - profile;
		
		if( opt.verbose && (mpi_rank == 0) ){
			cout << "Completed Bloom filter construction in " << profile << " sec" << endl;
		}
	}
	catch(const char *error){

        	cerr << "[" << mpi_rank << "] Caught the error: " << error << endl;
		
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        	return EXIT_FAILURE;
	}
	catch(...){

		cerr << "[" << mpi_rank << "] Caught an unhandled error" << endl;

		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		return EXIT_FAILURE;
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}

// Report our rank, the signal we caught and the time spent running
void terminate_program(int m_sig)
{
        cerr << "[" << mpi_rank << "] caught signal " << m_sig << endl;
        cerr << report_run_time() << endl;
        
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
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

