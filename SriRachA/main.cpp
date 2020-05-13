// Stream SRA sequence data
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <deque>
#include <unordered_map>
#include <algorithm>

#include <mpi.h>
#include <stdlib.h>
#include <getopt.h>
#include <signal.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif //_OPENMP

#include "sriracha.h"
#include "options.h"
#include "parse_sequence.h"
#include "word.h"
#include "sort.h"
#include "sra_stream.h"

using namespace std;

// MPI message tags
#define		RESULTS_LEN		1000
#define		RESULTS			1001

// Global variables for MPI
int mpi_numtasks;
int mpi_rank;

double start_time;

void terminate_program(int m_sig);
string report_run_time();

SRADownloadStatus search(const string &m_accession, ostream &m_out,
	const deque< pair< string, deque<Word> > > &m_subject_kmers, 
	StreamStats &m_info, const SrirachaOptions &m_opt);

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
		
		SrirachaOptions opt;
		ofstream fout;

		if(mpi_rank == 0){

			opt.load(argc, argv);

			// Currently only supporting search-by-kmer
			if( !opt.quit && (opt.search_strategy != SEARCH_BY_KMER) ){

				cerr << "Currently, SriRachA only supports search by kmer" << endl;
				opt.quit = true;
			}

			if( !opt.quit && !opt.output_filename.empty() ){
				
				fout.open( opt.output_filename.c_str() );

				if(!fout){

					cerr << "Unable to open " << opt.output_filename << " for writing" << endl;
					opt.quit = true;
				}
			}
		}

		broadcast(opt, mpi_rank, 0);

		if(opt.quit){

			MPI_Finalize();
			return EXIT_SUCCESS;
		}

		// All ranks other than 0 are silent
		if(mpi_rank > 0){
			opt.verbose = SILENT;
		}

		if(opt.verbose > TACITERN){
			cerr << "Running with " << mpi_numtasks << " MPI ranks" << endl;
		}

		ostream &out = ( (mpi_rank == 0) && fout.is_open() ) ? fout : cout;

		// Every rank loads the input (subject) sequences that will be searched
		// against each read (query) for every input accession.

		// For the search-by-kmer strategy
		deque< pair< string /*defline*/, deque<Word> /*kmers*/> > subject_kmers;

		for(deque<string>::const_iterator i = opt.input_sequence_files.begin();i != opt.input_sequence_files.end();++i){

			if(opt.verbose >= NORMAL){
				cerr << "Reading sequences from " << *i << endl;
			}

			SequenceIterator iter(*i);

			if(!iter){

				cerr << "[" << mpi_rank << "] Unable to open sequence file: " << *i << endl;
				throw __FILE__ ":main: Unable to open input sequence file";
			}

			while(iter){

				subject_kmers.push_back( make_pair( iter.get_info(), deque<Word>() ) );

				deque<Word> &local_kmers = subject_kmers.back().second;

				ForEachDuplexWord(iter.get_seq(), opt.kmer_len)

					if(ValidWord){
						local_kmers.push_back(CanonicalWord);
					}

				EndWord

				// Turn the list of extracted kmers into a set
				SORT( local_kmers.begin(), local_kmers.end() );
				local_kmers.erase( unique( local_kmers.begin(), local_kmers.end() ), local_kmers.end() );

				if(opt.verbose >= CHATTY){
					cerr << "\t" << iter.get_info() << " has " << local_kmers.size() << " unique kmers" << endl;
				}

				if( local_kmers.empty() ){

					if(opt.verbose >= TACITERN){
						cerr << "Did not extract any kmers from: " << iter.get_info() << endl;
					}

					subject_kmers.pop_back();
				}

				++iter;
			}
		}

		// Track if we have encounted an unrecoverable error
		bool failed_download = false;

		for(deque<string>::const_iterator a = opt.sra_accession.begin();a != opt.sra_accession.end();++a){

			if(opt.verbose >= NORMAL){
				cerr << "Searching " << *a << " ... ";
			}

			double download_profile = MPI_Wtime();

			StreamStats info;

			// Capture the SRADownloadStatus value so all MPI ranks can make sure that 
			// the entire SRA accession was correctly downloaded
			const SRADownloadStatus ret = search(*a, out, subject_kmers, info, opt);

			if( (ret != SRADownloadSuccess) && (mpi_rank == 0) ){

				cerr << "Unable to download SRA accession: " << *a 
					<< " (" << SRADownloadErrorStr[ret] << ")" << endl;
				
				failed_download = true;

				continue;
			}

			download_profile = MPI_Wtime() - download_profile;

			if(opt.verbose >= NORMAL){
				cerr << "complete in " << download_profile << " sec; " 
					<< info.num_reads << " reads and " << info.num_bases << " bases; "
					<< info.num_bases/(max(1.0, download_profile)*1.0e6) << " Mbp/sec" << endl;
			}
		}

		if( !opt.sra_accession_filename.empty() ){

			ifstream faccession( opt.sra_accession_filename.c_str() );

			if(!faccession){
				throw __FILE__ ":main: Unable to open filename with SRA accession list";
			}

			string accession;

			while(faccession >> accession){

				if(opt.verbose >= NORMAL){
					cerr << "Searching " << accession << " ... ";
				}

				double download_profile = MPI_Wtime();

				StreamStats info;

				// Capture the SRADownloadStatus value so all MPI ranks can make sure that 
				// the entire SRA accession was correctly downloaded
				const SRADownloadStatus ret = search(accession, out, subject_kmers, 
					info, opt);

				if( (ret != SRADownloadSuccess) && (mpi_rank == 0) ){

					cerr << "Unable to download SRA accession (from file): " << accession 
						<< " (" << SRADownloadErrorStr[ret] << ")" << endl;
					
					failed_download = true;

					continue;
				}

				download_profile = MPI_Wtime() - download_profile;

				if(opt.verbose >= NORMAL){

					cerr << "complete in " << download_profile << " sec; " 
						<< info.num_reads << " reads and " << info.num_bases << " bases; "
						<< info.num_bases/(max(1.0, download_profile)*1.0e6) << " Mbp/sec" << endl;
				}
			}
		}

		// Pipe from stdin on rank 0 if the user has not specified 
		// any other source of accessions
		if( opt.sra_accession.empty() && opt.sra_accession_filename.empty() ){

			if(opt.verbose >= TACITERN){
				cerr << "Reading SRA accession from stdin (ctrl+d to stop)" << endl;
			}

			while(true){

				string accession;

				if(mpi_rank == 0){

					if( !(cin >> accession) ){
						accession.clear();
					}
				}

				broadcast(accession, mpi_rank, 0);

				if( accession.empty() ){
					break;
				}

				if(opt.verbose >= NORMAL){
					cerr << "Searching " << accession<< " ... ";
				}

				double download_profile = MPI_Wtime();

				StreamStats info;

				// Capture the SRADownloadStatus value as an "int" so all MPI ranks can make sure that 
				// the entire SRA accession was correctly downloaded
				const SRADownloadStatus ret = search(accession, out, subject_kmers, 
					info, opt);

				if( (ret != SRADownloadSuccess) && (mpi_rank == 0) ){

					cerr << "Unable to download SRA accession (from stdin): " << accession 
						<< " (" << SRADownloadErrorStr[ret] << ")" << endl;
					
					failed_download = true;

					continue;
				}

				download_profile = MPI_Wtime() - download_profile;

				if(opt.verbose >= NORMAL){
					cerr << "complete in " << download_profile << " sec; " 
						<< info.num_reads << " reads and " << info.num_bases << " bases; "
						<< info.num_bases/(max(1.0, download_profile)*1.0e6) << " Mbp/sec" << endl;
				}
			}
		}

		const double profile = MPI_Wtime() - start_time;
		
		if(mpi_rank == 0){

			if( !failed_download ){

				// Indicate when *all* SRA accessions have been sucessfully searched.
				out << "//" << endl;
			}

			cout << "Completed SRA streaming in " << profile << " sec" << endl;
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


SRADownloadStatus search(const string &m_accession, ostream &m_out,
	const deque< pair< string, deque<Word> > > &m_subject_kmers, 
	StreamStats &m_info, const SrirachaOptions &m_opt)
{
	const size_t num_subject = m_subject_kmers.size();

	vector< deque<SearchMatch> > results(num_subject);
	vector<size_t> num_perfect_match(num_subject);

	// Store the SRA return codes as an "int" so that all ranks can
	// share their values with MPI_Reduce.
	int ret = SRADownloadSuccess;
	unsigned int attempt = 0;

	// Each rank will attempt to download the data and retry in the event of
	// a SRADownloadNetworkFailure.
	do{
		if(attempt > 0){

			// Reset the match results after every failed attempt
			results = vector< deque<SearchMatch> >(num_subject);
			num_perfect_match = vector<size_t>(num_subject);

			m_info = StreamStats();
		}

		switch(m_opt.search_strategy){
			case SEARCH_BY_ALIGN:
				throw __FILE__ ":search: SEARCH_BY_ALIGN not implemented yet";
				break;
			case SEARCH_BY_KMER:
				{
					void* param[] = {
						(void*)&results, 
						(void*)&m_subject_kmers, 
						(void*)&num_perfect_match,
						(void*)&m_opt, 
						(void*)NULL
					};

					ret = sra_stream(m_accession, search_by_kmer, param, &m_info);
				}
				break;
			case SEARCH_BY_BLOOM:
				throw __FILE__ ":search: SEARCH_BY_BLOOM not implemented yet";
				break;
			default:
				throw __FILE__ ":search: Invalid search strategy";
		};

		++attempt;
	}
	while( (ret == SRADownloadNetworkFailure) && (m_opt.max_retry > attempt) );

	// If a rank fails to download its portion of the accession, it will have
	// a ret value that is not equal to SRADownloadSuccess -- but don't act on this
	// yet! 

	// Cull the matches to avoid running out of RAM
	for(vector< deque<SearchMatch> >::iterator i = results.begin();i != results.end();++i){
			
		SORT( i->begin(), i->end() );

		if( (m_opt.max_num_match > 0) && (i->size() > m_opt.max_num_match) ){
			i->resize(m_opt.max_num_match);
		}
	}
	
	// Collect all of the search results for this accession from all ranks
	if(mpi_rank > 0){

		// Send the match results to the master
		unsigned int buffer_len = mpi_size(results);
		unsigned char* buffer = new unsigned char [buffer_len];

		if(buffer == NULL){
			throw __FILE__ ":search: Unable to allocate send buffer for results";
		}

		mpi_pack(buffer, results);

		if(MPI_Send( (void*)&buffer_len, 1, MPI_UNSIGNED, 0, RESULTS_LEN, MPI_COMM_WORLD ) != MPI_SUCCESS){
			throw __FILE__ ":search: Unable to send results buffer length";
		}

		if(MPI_Send( (void*)buffer, buffer_len, MPI_BYTE, 0, RESULTS, MPI_COMM_WORLD ) != MPI_SUCCESS){
			throw __FILE__ ":search: Unable to send results buffer";
		}

		delete [] buffer;
		buffer = NULL;
	}
	else{ // mpi_rank == 0

		for(int worker = 1;worker < mpi_numtasks;++worker){

			MPI_Status status;

			unsigned int buffer_len = 0;

			if(MPI_Recv(&buffer_len, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, RESULTS_LEN, MPI_COMM_WORLD, &status) != MPI_SUCCESS){
				throw __FILE__ ":search: Unable to receive buffer length";
			}

			unsigned char* buffer = new unsigned char[buffer_len];
			
			if(buffer == NULL){
				throw __FILE__ ":search: Unable to allocate receive buffer";
			}

			if(MPI_Recv(buffer, buffer_len, MPI_BYTE, status.MPI_SOURCE, RESULTS, MPI_COMM_WORLD, &status) != MPI_SUCCESS){
				throw __FILE__ ":search: Unable to receive buffer";
			}

			vector< deque<SearchMatch> > local_results;
			
			mpi_unpack(buffer, local_results);

			// Merge the local results in the results vector for rank 0
			for(size_t i = 0;i < num_subject;++i){
				
				results[i].insert( results[i].end(), 
					local_results[i].begin(), local_results[i].end() );
			}

			delete [] buffer;
			buffer = NULL;

			// Cull the matches to avoid running out of RAM
			for(vector< deque<SearchMatch> >::iterator i = results.begin();i != results.end();++i){
					
				SORT( i->begin(), i->end() );

				if( (m_opt.max_num_match > 0) && (i->size() > m_opt.max_num_match) ){
					i->resize(m_opt.max_num_match);
				}
			}
		}
	}

	// Did all of the MPI rank successfully download the accession? We need to check *before*
	// any output is written to disk!
	MPI_Allreduce(MPI_IN_PLACE, &ret, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	
	if(ret != SRADownloadSuccess){

		if(mpi_rank == 0){
			
			// Indicate in the output file that we were unable to search this accession
			m_out << m_accession << "\tNA\t0\t" << SRADownloadErrorStr[ret] << endl;
		}

		return SRADownloadStatus(ret);
	}

	// Collect some SRA file statistics
	MPI_Allreduce(MPI_IN_PLACE, &m_info.num_reads, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &m_info.num_bases, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

	// Only rank 0 outputs the results
	if(mpi_rank == 0){

		for(size_t i = 0;i < num_subject;++i){

			const deque<SearchMatch> &ref = results[i];

			// The results have already been sorted and culled to
			// have no more than opt.max_num_match elements
			// (unless opt.max_num_match == 0)
			for(deque<SearchMatch>::const_iterator j = ref.begin();j != ref.end();++j){

				m_out << m_accession << '\t' << j->read_index;
				
				if(j->read_subindex > 0){
					m_out << '.' << j->read_subindex;
				}
				
				m_out << '\t' << j->score 
					<< '\t' << m_subject_kmers[i].first << endl;
			}
		}
	}

	return SRADownloadStatus(ret);
}