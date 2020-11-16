#include <fstream>
#include <sstream>
#include <algorithm> // For lower_bound()

#include <math.h> // For ceil()

#include "maestro.h"
#include "sra_accession.h"
#include "binary_io.h"
#include "bloom.h"
#include "file_util.h"

// For openReadCollection (needed to estimate the number of bases in a SRA file)
#include <ncbi-vdb/NGS.hpp>

using namespace std;

extern int mpi_rank;
extern int mpi_numtasks;

#define	BYTES_PER_NODE	(16*GB)

bool parse_accession_loc(vector< pair<SraAccession, size_t /*location*/> > &m_accession_loc, const MaestroOptions &m_opt);
bool read_sra_repository(string &m_path);
bool restore_status(const string &m_filename, unsigned char* m_status, const size_t &m_num_sra, size_t &m_database_index);
bool write_status(const string &m_filename, unsigned char* m_status, const size_t &m_num_sra, size_t &m_database_index);
void display_status(ostream &m_out, unsigned char *m_status, const size_t &m_num_sra, const size_t &m_database_index,
	const size_t &m_num_kmer, const size_t &m_num_bp, const double &m_elapsed_time);
bool restore_download(deque<SraAccession> &m_download, 
	const string &m_sra_dir, 
	const vector< pair<SraAccession, size_t /*location*/> > &m_accession_loc,
	unsigned char *m_status, const bool &m_verbose);
bool restore_bloom(unordered_map< BloomParam, deque<SraAccession> > &m_bloom, const string &m_bloom_dir, 
	const vector< pair<SraAccession, size_t /*location*/> > &m_accession_loc, unsigned char *m_status,
	const bool &m_verbose);
bool process_event(MPI_Status &m_mpi_status, deque<int> &m_workers, unsigned char *m_status, 
	const vector< pair<SraAccession, size_t /*location*/> > &m_accession_loc, 
	deque<SraAccession> &m_download, deque<SraAccession> &m_retry, 
	const unsigned int &m_num_retry,
	unordered_map< BloomParam, deque<SraAccession> > &m_bloom, 
	size_t &m_num_kmer, size_t &m_num_bp,
	bool m_verbose);
size_t get_accession_index(const SraAccession &m_accession, 
	const vector< pair<SraAccession, size_t /*location*/> > &m_accession_loc);
bool schedule_database_creation(unordered_map< BloomParam, deque<SraAccession> > &m_bloom, 
	size_t &m_database_index, deque<int> &m_workers, const MAP<size_t, size_t> &m_database_size,
	bool m_verbose);
bool schedule_bloom_filter_creation(deque<SraAccession> &m_download, 
	deque<int> &m_workers, const vector< pair<SraAccession, size_t /*location*/> > &m_accession_loc, 
	const string &m_metadata_file, bool m_verbose);
bool schedule_streaming_bloom_filter_creation(size_t &m_curr_sra, deque<SraAccession> &m_retry, 
	deque<int> &m_workers, unsigned char* m_status, 
	const vector< pair<SraAccession, size_t /*location*/> > &m_accession_loc, 
	const string &m_metadata_file, const unsigned int &m_num_retry, bool m_verbose);
bool schedule_download(size_t &m_curr_sra, deque<SraAccession> &m_retry, deque<int> &m_workers, unsigned char* m_status, 
	const vector< pair<SraAccession, size_t> > &m_accession_loc, const unsigned int &m_num_retry, bool m_verbose);
vector<string> get_rank_names();

void maestro_main(MaestroOptions &m_opt)
{
	bool init = true;

	// Update the status information every five minutes
	const double status_update_every = 5*60;

	// Update the user every 15 minutes
	const double info_update_every = 15*60;

	if(m_opt.verbose){

		// Display the assay options
		cerr << m_opt << endl;
	}

	// Get the name of each MPI rank
	const vector<string> rank_name = get_rank_names();

	if(m_opt.verbose){
		
		for(int i = 0;i < mpi_numtasks;++i){
			cerr << "Rank " << i << " -> " << rank_name[i] << endl;
		}

		cerr << endl;
	}

	// To keep the database filesizes reasonable, the number of filters per database file
	// is adjusted based on the Bloom filter size
	MAP<size_t /*log 2 bloom filter size*/, size_t /* filters per database file*/> num_filters_per_database_file;

	for(size_t log_len = m_opt.min_log_2_filter_len;log_len <= m_opt.max_log_2_filter_len;++log_len){

		// If possible, we would like to have MAX_NUM_FILTER_CHUNK Bloom filters per database file (but no more).
		// However, we are not allowed to exceed MAX_DATABASE_FILE_SIZE_IN_GB gigabytes per database file!
		// In the equation for num bloom, 
		const size_t bits_per_byte = 8;
		const size_t num_bloom = (MAX_DATABASE_FILE_SIZE_IN_GB * bits_per_byte * GB)/(1UL << log_len);

		num_filters_per_database_file[log_len] = min(MAX_NUM_FILTER_CHUNK, num_bloom);

		if(m_opt.verbose){

			cerr << "Number of L=" << log_len << " Bloom filters per database file = " 
				<< num_filters_per_database_file[log_len] << endl;
		}
	}

	// Step 1: Make sure that the user has specified a local repository for SRA files.
	// We need the location of this file repository to be able to clean up SRA files
	// during the course of database creation (to avoid running out of space).
	string local_sra_repository_dir;

	if( !read_sra_repository(local_sra_repository_dir) ){

		init = false;
		broadcast(init, mpi_rank, 0);
		return;
	}

	if(m_opt.verbose){
		cerr << "Local SRA file repository = " << local_sra_repository_dir << endl;
	}

	// Step 2: Make sure that the scratch directory exits
	// To make file cleanup easier, we will create (if they do not already exist) the following
	// subdirectories:
	//		scratch_dir/
	//			bloom/
	//			database/
	if(!is_dir(m_opt.scratch_bloom_dir)){
		if( !make_dir(m_opt.scratch_bloom_dir) ){

			cerr << "Unable to create the scratch bloom directory: " << m_opt.scratch_bloom_dir << endl;

			init = false;
			broadcast(init, mpi_rank, 0);
			return;
		}
	}

	if(!is_dir(m_opt.scratch_database_dir)){
		if( !make_dir(m_opt.scratch_database_dir) ){

			cerr << "Unable to create the scratch database directory:" << m_opt.scratch_database_dir << endl;

			init = false;
			broadcast(init, mpi_rank, 0);
			return;
		}
	}
	
	// Step 3:  Read the binary metadata file and extract:
	//	- The list of SRA accessions
	//	- The starting location of the each SRA accession metadata record in the metadata file
	// Note that the resulting list of (accession, location) will be sorted by accession for
	// fast lookup by accession
	vector< pair<SraAccession, size_t /*location*/> > accession_loc;

	init = parse_accession_loc(accession_loc, m_opt.metadata_file, m_opt.verbose);

	if(!init){

		broadcast(init, mpi_rank, 0);
		return;
	}

	const size_t num_sra = accession_loc.size();

	if(num_sra == 0){

		init = false;

		cerr << "Did not read any SRA accessions from the input metadata file" << endl;

		broadcast(init, mpi_rank, 0);
		return;
	}

	// Record the progress of the search in an array of num_sra bytes
	size_t database_index = 1; // <-- Start counting database files from 1
	unsigned char *status = new unsigned char [num_sra];

	if(!status){
		throw __FILE__ ":maestro_main: Unable to allocate status buffer";
	}

	memset(status, STATUS_INIT, num_sra);

	init = restore_status(m_opt.status_file, status, num_sra, 
		database_index, true /*create missing*/);

	if(!init){

		cerr << "Unable to restore from status file" << endl;

		broadcast(init, mpi_rank, 0);
		return;
	}

	if(m_opt.retry_bloom){

		// Force the inclusion of all Bloom filter failures (including the final failure state
		// of STATUS_BLOOM_FAIL).
		for(size_t i = 0;i < num_sra;++i){

			switch(status[i]){

				case STATUS_BLOOM_FAIL:
				case STATUS_BLOOM_FAIL_1: case STATUS_BLOOM_FAIL_2: case STATUS_BLOOM_FAIL_3:
				case STATUS_BLOOM_FAIL_4: case STATUS_BLOOM_FAIL_5: case STATUS_BLOOM_FAIL_6:
				case STATUS_BLOOM_FAIL_7: case STATUS_BLOOM_FAIL_8: case STATUS_BLOOM_FAIL_9:
				case STATUS_BLOOM_FAIL_10:
					status[i] = STATUS_INIT;
					break;
			};
		}
	}

	for(deque<SraAccession>::const_iterator i = m_opt.skip_sra.begin();i != m_opt.skip_sra.end();++i){

		try{

			const size_t index = get_accession_index(*i, accession_loc);

			status[index] = STATUS_SKIPPED;

			if(m_opt.verbose){
				cerr << "Skipping " << accession_to_str(*i) << endl;
			}
		}
		catch(...){

			cerr << "A user-specified SRA run accession to skip, " << 
				accession_to_str(*i) << ", is not valid" << endl;
			init = false;

			broadcast(init, mpi_rank, 0);
			return;
		}
	}

	if(m_opt.verbose){
		display_status(cerr, status, num_sra, database_index, 
			0 /*running kmer*/, 0 /*running bp*/, -1.0 /*elapsed sec*/);
	}

	// Queues to track the SRA records that are actively being processed
	deque<SraAccession> sra_retry; // Needs to be retried (either download or streaming Bloom)
	deque<SraAccession> sra_download; // Successfully downloaded, needs Bloom conversion
	unordered_map< BloomParam, deque<SraAccession> > sra_bloom; // Successfull Bloom filter, needs database

	// When directly streaming SRA data, we don't need to restore downloads
	if(m_opt.stream_sra == false){

		init = restore_download(sra_download, local_sra_repository_dir, accession_loc, status, 
			m_opt.verbose);

		if(!init){

			cerr << "Unable to restore downloads" << endl;

			broadcast(init, mpi_rank, 0);
			return;
		}
	}
	
	init = restore_bloom(sra_bloom, m_opt.scratch_bloom_dir, accession_loc, status, m_opt.verbose);

	if(!init){

		cerr << "Unable to restore Bloom filters" << endl;

		broadcast(init, mpi_rank, 0);
		return;
	}

	// Let the other ranks know that the initialization was successfull
	broadcast(init, mpi_rank, 0);

	broadcast(local_sra_repository_dir, mpi_rank, 0);
	
	size_t curr_sra = 0;
	size_t end_sra = num_sra;

	if(m_opt.limit_num_download > 0){

		// The number of SRA records that still need to be downloaded
		size_t count = 0;

		for(size_t i = curr_sra;i < num_sra;++i){

			switch(status[i]){
				case STATUS_INIT:
				case STATUS_DOWNLOAD_FAIL_1: case STATUS_DOWNLOAD_FAIL_2: case STATUS_DOWNLOAD_FAIL_3:
				case STATUS_DOWNLOAD_FAIL_4: case STATUS_DOWNLOAD_FAIL_5: case STATUS_DOWNLOAD_FAIL_6:
				case STATUS_DOWNLOAD_FAIL_7: case STATUS_DOWNLOAD_FAIL_8: case STATUS_DOWNLOAD_FAIL_9:
				case STATUS_DOWNLOAD_FAIL_10:
				//case STATUS_BLOOM_FAIL: <-- Don't include hard Bloom filter failures
				case STATUS_BLOOM_FAIL_1: case STATUS_BLOOM_FAIL_2: case STATUS_BLOOM_FAIL_3:
				case STATUS_BLOOM_FAIL_4: case STATUS_BLOOM_FAIL_5: case STATUS_BLOOM_FAIL_6:
				case STATUS_BLOOM_FAIL_7: case STATUS_BLOOM_FAIL_8: case STATUS_BLOOM_FAIL_9:
				case STATUS_BLOOM_FAIL_10:
					++count;
					break;
			};

			if(count == m_opt.limit_num_download){

				end_sra = i + 1;
				break;
			}
		}

		if(m_opt.verbose){
			cerr << "Limiting number of SRA downloads to " << count << " records" << endl;
		}
	}
	
	// The list of free workers
	deque<int> workers;

	for(int w = 1;w < mpi_numtasks;++w){
		workers.push_back(w);
	}

	// The main event loop
	MPI_Status mpi_status;
	int message_ready;
	double last_status_update = MPI_Wtime();
	double last_info_update = last_status_update;
	double last_download = last_status_update;
	size_t running_num_kmer = 0;
	size_t running_num_bp = 0;

	if(m_opt.verbose){
		cerr << "Entering main event loop with " << workers.size() << " workers" << endl;
	}

	while(true){

		// Are we done?
		if( (curr_sra >= end_sra) && sra_download.empty() && sra_bloom.empty() &&
			sra_retry.empty() && ( (int)( workers.size() ) == (mpi_numtasks - 1) ) ){

			cerr << "Database construction complete!" << endl;
			break;
		}

		// Update the status on disk if needed
		if( (MPI_Wtime() - last_status_update) > status_update_every){

			if( !write_status(m_opt.status_file, status, num_sra, database_index) ){
				cerr << "** Warning! Unable to commit status information to disk!" << endl;
			}

			last_status_update = MPI_Wtime();
		}

		if( m_opt.verbose && (MPI_Wtime() - last_info_update > info_update_every) ){

			const size_t max_num_workers = mpi_numtasks - 1;

			cerr << workers.size() << "/" << max_num_workers << " workers are idle" << endl;
			cerr << (max_num_workers - workers.size()) << "/" << max_num_workers
				<< " workers are actively processing" << endl;
				
			display_status(cerr, status, num_sra, database_index, running_num_kmer, running_num_bp,
				MPI_Wtime() - last_info_update);

			running_num_kmer = 0;
			running_num_bp = 0;

			last_info_update = MPI_Wtime();
		}

		message_ready = 0;

		// Check for incomming messages
		if(MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &message_ready, &mpi_status) != MPI_SUCCESS){
			throw __FILE__ ":maestro_main: Error in MPI_Iprobe";
		}

		// Handle incomming messages
		if(message_ready){

			const bool save_status = process_event(mpi_status, workers, status, accession_loc, 
				sra_download, sra_retry, m_opt.num_download_attempt, 
				sra_bloom, running_num_kmer, running_num_bp, m_opt.verbose);
			
			// If we successfull sent a database file to S3, then update the status on disk
			if(save_status){

				if( !write_status(m_opt.status_file, status, num_sra, database_index) ){
					cerr << "** Warning! Unable to commit the status information to disk after process_even!" << endl;
				}

				last_status_update = MPI_Wtime();
			}

			continue;
		}

		// Priority #1: Creating database files and sending them to S3 storage. For now, require that
		// there are filters_per_database_file Bloom filters that all share the same Bloom parameters before
		// building a database. This should probably a user parameter.
		//
		// If there are no more SRA records to download and no more Bloom filters to create *and* all of the
		// workers are idle, then we need to pack the remaining Bloom filters into database files (even if there are
		// fewer than filters_per_database_file Bloom filters per Bloom parameter set)
		if( (curr_sra >= end_sra) && sra_download.empty() && ( (int)( workers.size() ) == (mpi_numtasks - 1) ) ){

			for(size_t log_len = m_opt.min_log_2_filter_len;log_len <= m_opt.max_log_2_filter_len;++log_len){
				num_filters_per_database_file[log_len] = 0; // <-- zero is a special value that *forces* database file construction
			}
		}

		if( schedule_database_creation(sra_bloom, database_index, workers, num_filters_per_database_file, m_opt.verbose) ){
			continue;
		}

		if(m_opt.stream_sra){

			if( (MPI_Wtime() - last_download) < m_opt.download_delay ){
				continue;
			}

			// Priority #2: Creating Bloom filters for SRA records by directly streaming the SRA data
			if( schedule_streaming_bloom_filter_creation(curr_sra, sra_retry, workers, status, accession_loc, 
				m_opt.metadata_file, m_opt.num_download_attempt, m_opt.verbose) ){

				last_download = MPI_Wtime();
				continue;
			}
		}
		else{

			// Priority #2: Creating Bloom filters for SRA records that have already been downloaded
			if( schedule_bloom_filter_creation(sra_download, workers, accession_loc, m_opt.metadata_file, m_opt.verbose) ){
				continue;
			}

			if( (MPI_Wtime() - last_download) < m_opt.download_delay ){
				continue;
			}

			// Priority #3: Downloading new SRA records to the local repository. Only initiate a download if we have
			// not exceeded the requested number of SRA records to download (i.e. the end_sra value that is either num_sra
			// or a limit imposed by the user in m_opt.limit_num_download).
			if( (curr_sra < end_sra) && 
				schedule_download(curr_sra, sra_retry, workers, status, accession_loc, m_opt.num_download_attempt, 
					m_opt.verbose) ){
				
				last_download = MPI_Wtime();
				continue;
			}
		}
	}

	// Make a final update of the status information on disk
	if( !write_status(m_opt.status_file, status, num_sra, database_index) ){
		cerr << "** Warning! Unable to commit the last status information to disk!" << endl;
	}

	if( m_opt.verbose){
		
		// Final status report
		display_status(cerr, status, num_sra, database_index, running_num_kmer, running_num_bp,
			MPI_Wtime() - last_info_update);
	}

	// Tell all of the workers to shut down
	for(int w = 1;w < mpi_numtasks;++w){

		// Send a zero byte message to tell the workers to quit
		if( MPI_Send( NULL, 0, MPI_BYTE, w, MAESTRO_QUIT, MPI_COMM_WORLD ) != MPI_SUCCESS ){
			throw __FILE__ ":maestro_main: Error sending MAESTRO_QUIT to worker";
		}
	}

	if(status){
		delete [] status;
	}
}

void display_status(ostream &m_out, unsigned char *m_status, const size_t &m_num_sra, const size_t &m_database_index,
	const size_t &m_num_kmer, const size_t &m_num_bp, const double &m_elapsed_time){

	// Summarize the current state of the SRA database creation process
	unordered_map<unsigned char, size_t> status_count;

	for(size_t i = 0;i < m_num_sra;++i){
		++status_count[ m_status[i] ];
	}

	const time_t current_time = time(NULL);
	const char* time_stamp = ctime(&current_time);

	m_out << "Current database index is " << m_database_index << endl;

	if(time_stamp != NULL){
		m_out << "Current SRA database creation status @ " << time_stamp;
	}
	else{
		m_out << "Current SRA database creation status:" << endl;
	}

	if(m_elapsed_time > 0.0){
		m_out << "Bloom filter rate is " << m_num_kmer/m_elapsed_time << " kmers/sec and " 
			<< m_num_bp/m_elapsed_time << " bp/sec" << endl;
	}

	if(status_count[STATUS_INIT] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_INIT])/m_num_sra 
			<< "% (" << status_count[STATUS_INIT] << ") records need to be processed" << endl;
	}

	if(status_count[STATUS_DOWNLOAD_SUCCESS] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DOWNLOAD_SUCCESS])/m_num_sra 
			<< "% (" << status_count[STATUS_DOWNLOAD_SUCCESS] << ") records have been successfully downloaded" << endl;
	}

	if(status_count[STATUS_DOWNLOAD_FAIL] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DOWNLOAD_FAIL])/m_num_sra 
			<< "% (" << status_count[STATUS_DOWNLOAD_FAIL] << ") records could not be downloaded" << endl;
	}

	if(status_count[STATUS_DOWNLOAD_FAIL_1] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DOWNLOAD_FAIL_1])/m_num_sra 
			<< "% (" << status_count[STATUS_DOWNLOAD_FAIL_1] << ") records could not be downloaded after 1 attempt" << endl;
	}

	if(status_count[STATUS_DOWNLOAD_FAIL_2] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DOWNLOAD_FAIL_2])/m_num_sra 
			<< "% (" << status_count[STATUS_DOWNLOAD_FAIL_2] << ") records could not be downloaded after 2 attempts" << endl;
	}

	if(status_count[STATUS_DOWNLOAD_FAIL_3] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DOWNLOAD_FAIL_3])/m_num_sra 
			<< "% (" << status_count[STATUS_DOWNLOAD_FAIL_3] << ") records could not be downloaded after 3 attempts" << endl;
	}

	if(status_count[STATUS_DOWNLOAD_FAIL_4] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DOWNLOAD_FAIL_4])/m_num_sra 
			<< "% (" << status_count[STATUS_DOWNLOAD_FAIL_4] << ") records could not be downloaded after 4 attempts" << endl;
	}

	if(status_count[STATUS_DOWNLOAD_FAIL_5] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DOWNLOAD_FAIL_5])/m_num_sra 
			<< "% (" << status_count[STATUS_DOWNLOAD_FAIL_5] << ") records could not be downloaded after 5 attempts" << endl;
	}

	if(status_count[STATUS_DOWNLOAD_FAIL_6] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DOWNLOAD_FAIL_6])/m_num_sra 
			<< "% (" << status_count[STATUS_DOWNLOAD_FAIL_6] << ") records could not be downloaded after 6 attempts" << endl;
	}

	if(status_count[STATUS_DOWNLOAD_FAIL_7] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DOWNLOAD_FAIL_7])/m_num_sra 
			<< "% (" << status_count[STATUS_DOWNLOAD_FAIL_7] << ") records could not be downloaded after 7 attempts" << endl;
	}

	if(status_count[STATUS_DOWNLOAD_FAIL_8] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DOWNLOAD_FAIL_8])/m_num_sra 
			<< "% (" << status_count[STATUS_DOWNLOAD_FAIL_8] << ") records could not be downloaded after 8 attempts" << endl;
	}

	if(status_count[STATUS_DOWNLOAD_FAIL_9] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DOWNLOAD_FAIL_9])/m_num_sra 
			<< "% (" << status_count[STATUS_DOWNLOAD_FAIL_9] << ") records could not be downloaded after 9 attempts" << endl;
	}

	if(status_count[STATUS_DOWNLOAD_FAIL_10] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DOWNLOAD_FAIL_10])/m_num_sra 
			<< "% (" << status_count[STATUS_DOWNLOAD_FAIL_10] << ") records could not be downloaded after 10 attempts" << endl;
	}

	if(status_count[STATUS_DOWNLOAD_FAIL] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DOWNLOAD_FAIL])/m_num_sra 
			<< "% (" << status_count[STATUS_DOWNLOAD_FAIL] << ") records could not be downloaded (gave up)" << endl;
	}

	if(status_count[STATUS_BLOOM_SUCCESS] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_BLOOM_SUCCESS])/m_num_sra 
			<< "% (" << status_count[STATUS_BLOOM_SUCCESS] << ") records have been stored in Bloom filters" << endl;
	}

	if(status_count[STATUS_BLOOM_INVALID] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_BLOOM_INVALID])/m_num_sra 
			<< "% (" << status_count[STATUS_BLOOM_INVALID] << ") records could not be converted into Bloom filters (FP rate or no data)" << endl;
	}

	if(status_count[STATUS_BLOOM_FAIL] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_BLOOM_FAIL])/m_num_sra 
			<< "% (" << status_count[STATUS_BLOOM_FAIL] << ") records could not be stored in Bloom filters (error)" << endl;
	}

	if(status_count[STATUS_DATABASE_SUCCESS] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DATABASE_SUCCESS])/m_num_sra 
			<< "% (" << status_count[STATUS_DATABASE_SUCCESS] << ") records were stored in database files" << endl;
	}

	if(status_count[STATUS_DATABASE_FAIL] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DATABASE_FAIL])/m_num_sra 
			<< "% (" << status_count[STATUS_DATABASE_FAIL] << ") records could not be stored in database files" << endl;
	}

	if(status_count[STATUS_DATABASE_UPLOAD_FAIL] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_DATABASE_UPLOAD_FAIL])/m_num_sra 
			<< "% (" << status_count[STATUS_DATABASE_UPLOAD_FAIL] << ") records could not be uploaded to S3" << endl;
	}

	if(status_count[STATUS_SKIPPED] > 0){
		m_out << "\t" << (100.0*status_count[STATUS_SKIPPED])/m_num_sra 
			<< "% (" << status_count[STATUS_SKIPPED] << ") records were skipped" << endl;
	}
}

bool restore_download(deque<SraAccession> &m_download, 
	const string &m_sra_dir, 
	const vector< pair<SraAccession, size_t /*location*/> > &m_accession_loc,
	unsigned char *m_status, const bool &m_verbose)
{
	// Always return true 
	const bool ret = true;

	size_t num_restore = 0;
	const size_t num_sra = m_accession_loc.size();

	for(size_t i = 0;i < num_sra;++i){

		if(m_status[i] == STATUS_DOWNLOAD_SUCCESS){

			const string acc = accession_to_str(m_accession_loc[i].first);

			// Make sure we have a valid .sra file in the download directory
			if( ! is_file(m_sra_dir + PATH_SEPARATOR + acc + ".sra") ){
				
				cerr << "Unable to find SRA file for accession " << acc << endl;

				// Since we couldn't find the SRA data, roll back the status
				// to STATUS_INIT and try again.
				m_status[i] = STATUS_INIT;

				continue;
			}

			m_download.push_back(m_accession_loc[i].first);
		}
	}

	if( ret && m_verbose && (num_restore > 0) ){
		cerr << "Restored a total of " << num_restore << " SRA downloads" << endl;
	}

	return ret;
}

size_t estimate_num_bases(const string &m_accession)
{
	// The openReadCollection() will first look for a local SRA dataset for the
	// provided accession. If this fails, it will attempt to download. Since we
	// are using a local repository for SRA files, this should always find a local
	// repository.
	ngs::ReadCollection run(  ncbi::NGS::openReadCollection(m_accession) );
		
	// Note that num_read is the number of either paired or
	// unpaired reads. For paired reads, this is half the
	// the number of sequences!
	//
	// We will estimate the number of bases by reading either
	// run.getReadCount(ngs::Read::all) or SAMPLE_READS, whichever
	// is smaller
	#define		SAMPLE_READS		10000

	const size_t num_total_read = run.getReadCount(ngs::Read::all);
	const size_t num_sample_read = min(size_t(SAMPLE_READS), num_total_read);
						
	ngs::ReadIterator run_iter = 
		ngs::ReadIterator( run.getReadRange ( 1, num_sample_read, ngs::Read::all ) );
		
	// Start read counting at zero, since we will be adding this to "start"
	size_t base_count = 0;
					
	while( run_iter.nextRead() ){
			
		while( run_iter.nextFragment() ){
			base_count += run_iter.getFragmentBases().toString().size();
		}
	}

	// Use the sampled reads to estimate the total number of bases
	if(num_sample_read > 0){
		base_count = size_t( base_count*(double(num_total_read)/num_sample_read) );
	}

	return base_count;
}

bool restore_bloom(unordered_map< BloomParam, deque<SraAccession> > &m_bloom, const string &m_bloom_dir, 
	const vector< pair<SraAccession, size_t /*location*/> > &m_accession_loc, unsigned char *m_status, 
	const bool &m_verbose)
{
	// If we encounter an error for a given Bloom filter, we will modify the status
	// to ensure that the Bloom filter is recomputed
	const bool ret = true; // Always return true (for now).
	size_t num_restore = 0;
	size_t num_error = 0;

	const size_t num_sra = m_accession_loc.size();

	for(size_t i = 0;i < num_sra;++i){

		// Recover both types of vaild Bloom filter files
		if( (m_status[i] == STATUS_BLOOM_SUCCESS) || 
		    (m_status[i] == STATUS_DATABASE_FAIL) ){ // <-- Don't try to recover STATUS_DATABASE_UPLOAD_FAIL
													 // since the Bloom filters will be gone! These files
													 // must be manually uploaded to the S3 bucket.

			const string acc = accession_to_str(m_accession_loc[i].first);
			const string filename = m_bloom_dir + PATH_SEPARATOR + acc + ".bloom";

			// Make sure we have a valid .bloom file in the download directory
			if( ! is_file(filename) ){
				
				cerr << "Unable to find Bloom filter file for accession " << acc << endl;
				++num_error;

				// Roll back the status to download the SRA and recompute this Bloom filter
				m_status[i] = STATUS_INIT;
				continue;
			}

			ifstream fin(filename.c_str(), ios::binary);

			if(!fin){

				cerr << "Unable to open Bloom filter " << filename << " for reading" << endl;
				
				++num_error;

				// Roll back the status to download the SRA and recompute this Bloom filter
				m_status[i] = STATUS_INIT;

				continue;
			}

			// Read the Bloom filter parameters from the Bloom filter file header
			try{
				
				unsigned char magic;

				binary_read(fin, magic);

				if(magic == BLOOM_MAGIC_COMPLETE){

					BloomParam param;

					binary_read(fin, param);

					m_bloom[param].push_back(m_accession_loc[i].first);

					if(m_verbose){
						cerr << "Restored Bloom filter for " << acc << endl;	
					}

					++num_restore;
				}
				else{

					cerr << "The completion bit is not set for " << acc << endl;
					
					++num_error;

					// Roll back the status to download the SRA and recompute this Bloom filter
					m_status[i] = STATUS_INIT;
				}
			}
			catch(...){

				cerr << "Error trying to read the header for " << acc << endl;

				++num_error;

				// Roll back the status to download the SRA and recompute this Bloom filter
				m_status[i] = STATUS_INIT;
			}
		}
	}

	if(m_verbose){

		if( ret && (num_restore > 0) ){
			cerr << "Restored a total of " << num_restore << " Bloom filters" << endl;
		}

		if(num_error > 0){
			cerr << "Unable to restore " << num_error << " Bloom filters. These will be recomputed" << endl;
		}
	}

	return ret;
}

// Return true if we need to immediately update the status on disk
// Return false otherwise
bool process_event(MPI_Status &m_mpi_status, deque<int> &m_workers, 
	unsigned char *m_status, 
	const vector< pair<SraAccession, size_t /*location*/> > &m_accession_loc, 
	deque<SraAccession> &m_download, deque<SraAccession> &m_retry, 
	const unsigned int &m_num_retry,
	unordered_map< BloomParam, deque<SraAccession> > &m_bloom,
	size_t &m_num_kmer, size_t &m_num_bp,
	bool m_verbose)
{
	bool ret = false;

	int buffer_size = 0;
	
	// All messages are packed into byte arrays
	if(MPI_Get_count(&m_mpi_status, MPI_BYTE, &buffer_size) != MPI_SUCCESS){
		throw __FILE__ ":process_event: Error obtaining main even buffer size with MPI_Get_count";
	}

	if(buffer_size <= 0){
		throw __FILE__ ":process_event: Invalid main even buffer size";
	}

	unsigned char *buffer = new unsigned char [buffer_size];

	if(buffer == NULL){
		throw __FILE__ ":process_event: Unable to allocate main even buffer";
	}

	if(MPI_Recv(buffer, buffer_size, MPI_BYTE, m_mpi_status.MPI_SOURCE, 
		m_mpi_status.MPI_TAG, MPI_COMM_WORLD, &m_mpi_status) != MPI_SUCCESS){

		throw __FILE__ ":process_event: Unable to receive main even buffer";
	}

	unsigned char* ptr = buffer;

	SraAccession acc;
	vector<SraAccession> acc_list;
	BloomParam param;
	BloomProgress progress;
	size_t index;
	double profile_sec;
	int prefetch_error;
	float worker_mem;

	switch(m_mpi_status.MPI_TAG){

		case STATUS_DOWNLOAD_SUCCESS:
			
			// Receive the:
			//	- accession 
			//  - time spent downloading
			//	- %worker memory used
			ptr = mpi_unpack(ptr, acc);
			ptr = mpi_unpack(ptr, profile_sec);
			ptr = mpi_unpack(ptr, worker_mem);
			
			m_download.push_back(acc);

			// Update the accession status in RAM
			m_status[get_accession_index(acc, m_accession_loc)] = m_mpi_status.MPI_TAG;

			// Return the worker to the pool of available workers
			m_workers.push_back(m_mpi_status.MPI_SOURCE);

			if(m_verbose){

				cerr << "Worker " << m_mpi_status.MPI_SOURCE << " successfully downloaded " 
					<< accession_to_str(acc) << " in " << profile_sec << " sec; mem = " 
					<< worker_mem << "%" << endl;
			}

			break;

		case STATUS_DOWNLOAD_FAIL:

			// Receive the:
			//	- accession
			//  - time spent downloading
			//  - prefetch error code
			//	- %worker memory used
			ptr = mpi_unpack(ptr, acc);
			ptr = mpi_unpack(ptr, profile_sec);
			ptr = mpi_unpack(ptr, prefetch_error);
			ptr = mpi_unpack(ptr, worker_mem);

			// Increment the failure counter (which is stored in the status array)
			index = get_accession_index(acc, m_accession_loc);

			if(m_status[index] == STATUS_INIT){
				m_status[index] = STATUS_DOWNLOAD_FAIL_1;
			}
			else{
				++m_status[index];
			}

			if(m_status[index] > (STATUS_DOWNLOAD_FAIL + m_num_retry) ){
				m_status[index] = STATUS_DOWNLOAD_FAIL;
			}
			else{

				// Put this accession in the retry queue
				m_retry.push_back(acc);
			}

			// Return the worker to the pool of available workers
			m_workers.push_back(m_mpi_status.MPI_SOURCE);

			if(m_verbose){

				cerr << "Worker " << m_mpi_status.MPI_SOURCE << " failed to download " 
					<< accession_to_str(acc) << " in " << profile_sec << " sec; mem = " 
					<< worker_mem << "%;";
				
				if(m_status[index] == STATUS_DOWNLOAD_FAIL){
					cerr << "final";
				}
				else{
					cerr << "attempt " << int(m_status[index] - STATUS_DOWNLOAD_FAIL_1) - 1;
				}

				cerr << "; prefetch error 0x" 
					<< std::hex << prefetch_error << std::dec << endl;
			}

			break;

		case STATUS_BLOOM_SUCCESS:

			// Receive the:
			//	- accession
			//	- Bloom param
			//	- Bloom progress
			//	- time to build Bloom filter in sec
			//	- %worker memory used
			ptr = mpi_unpack(ptr, acc);
			ptr = mpi_unpack(ptr, param);
			ptr = mpi_unpack(ptr, progress);
			ptr = mpi_unpack(ptr, profile_sec);
			ptr = mpi_unpack(ptr, worker_mem);

			// Schedule this Bloom filter for database inclusion
			m_bloom[param].push_back(acc);

			// Update the accession status in RAM
			m_status[get_accession_index(acc, m_accession_loc)] = m_mpi_status.MPI_TAG;

			// Return the worker back to the active pool
			m_workers.push_back(m_mpi_status.MPI_SOURCE);

			// Track the number of kmers and the number of bases that have
			// been stored in Bloom filters
			m_num_kmer += progress.num_kmer;
			m_num_bp += progress.num_bp;

			if(m_verbose){
				
				// deflation == (number of Bloom filter bytes)/(sequence bytes)
				// We assume that:
				//	- The number of bases in an SRA record == num_bp
				//	- Each base requires two bits (4 bases per byte)
				//
				// deflation = (number of filter bits/8)/(num_bp/4)
				//           = (number of filter bits)/(2*num_bp)
				const float deflation = (progress.num_bp == 0) ? 0.0f:
					float( param.filter_len() )/(2*progress.num_bp);

				// Track the ratio of kmers/bp -- how much repeated sequence is there in
				// each SRA record? Restated: what fraction of the kmers are unique?
				const float uniqueness = (progress.num_bp == 0) ? 0.0f:
					float(progress.num_kmer)/progress.num_bp;

				cerr << "Worker " << m_mpi_status.MPI_SOURCE << " created Bloom filter for " << accession_to_str(acc) 
					<< " in " << profile_sec << " sec; mem = " << worker_mem 
					<< "%; L = " << param.log_2_filter_len 
					<< "; count L = " << progress.log_2_counting_filter_len
					<< "; h = " << param.num_hash << "; unique = " << uniqueness 
					<< "; deflation = " << deflation << endl;
			}

			break;

		case STATUS_BLOOM_FAIL:

			// Receive the:
			//	- accession
			//	- Bloom progress
			//	- time spend attempting to make Bloom filter
			//	- %worker memory used
			ptr = mpi_unpack(ptr, acc);
			ptr = mpi_unpack(ptr, progress);
			ptr = mpi_unpack(ptr, profile_sec);
			ptr = mpi_unpack(ptr, worker_mem);

			// Increment the failure counter (which is stored in the status array)
			index = get_accession_index(acc, m_accession_loc);

			// Update the accession status in RAM
			if( (m_status[index] == STATUS_INIT) || (m_status[index] == STATUS_DOWNLOAD_SUCCESS) ){
				m_status[index] = STATUS_BLOOM_FAIL_1;
			}
			else{
				++m_status[index];
			}

			// Note that STATUS_BLOOM_FAIL_1 != STATUS_BLOOM_FAIL, hence the need for
			// the following awkward test.
			if( (m_status[index] + 1U) > (unsigned int)(STATUS_BLOOM_FAIL_1 + m_num_retry) ){

				// A "hard" failure -- stop retrying!
				m_status[index] = STATUS_BLOOM_FAIL;
			}
			else{

				// Put this accession in the retry queue
				m_retry.push_back(acc);
			}

			// Return the worker back to the active pool
			m_workers.push_back(m_mpi_status.MPI_SOURCE);

			if(m_verbose){

				cerr << "Worker " << m_mpi_status.MPI_SOURCE << " failed to create Bloom filter for " 
					<< accession_to_str(acc) << " in " << profile_sec << " sec; mem = " << worker_mem << "%;";
				
				if(m_status[index] == STATUS_BLOOM_FAIL){
					cerr << "final";
				}
				else{
					cerr << "attempt " << int(m_status[index] - STATUS_BLOOM_FAIL_1) - 1;
				}

				cerr << endl;

				// DEBUG -- how far does streaming go before failing?

				cerr << "\tvalid_read_collection = " << (progress.valid_read_collection ? "true" : "false") << endl;

				if(progress.num_primary_align > 0){

					cerr << "\tPrimary alignment: " << progress.curr_primary_align << " out of " 
						<< progress.num_primary_align << " (" 
						<< (100.0*progress.curr_primary_align)/progress.num_primary_align << "%)" << endl;
				}

				if(progress.num_unaligned_read > 0){

					cerr << "\tUnaligned reads: " << progress.curr_unaligned_read << " out of " 
						<< progress.num_unaligned_read << " (" 
						<< (100.0*progress.curr_unaligned_read)/progress.num_unaligned_read << "%)" << endl;

					cerr << "\tCurrent fragment: " << progress.curr_fragment << endl;
				}

				if(progress.num_read > 0){

					cerr << "\tReads: " << progress.curr_read << " out of " 
						<< progress.num_read << " (" << (100.0*progress.curr_read)/progress.num_read << "%)" << endl;
					cerr << "\tCurrent fragment: " << progress.curr_fragment << endl;
				}

				if( !progress.error.empty() ){
					cerr << "\tError: " << progress.error << endl;
				}
			}

			break;

		case STATUS_BLOOM_INVALID:

			// Receive the:
			//	- accession
			//	- Bloom progress
			//	- time spend attempting to build Bloom filter
			//	- %worker memory used
			ptr = mpi_unpack(ptr, acc);
			ptr = mpi_unpack(ptr, progress);
			ptr = mpi_unpack(ptr, profile_sec);
			ptr = mpi_unpack(ptr, worker_mem);

			// Update the accession status in RAM
			m_status[get_accession_index(acc, m_accession_loc)] = m_mpi_status.MPI_TAG;

			// Return the worker back to the active pool
			m_workers.push_back(m_mpi_status.MPI_SOURCE);

			if(m_verbose){

				cerr << "Worker " << m_mpi_status.MPI_SOURCE << " cannot find valid Bloom filter parameters for " 
					<< accession_to_str(acc) << " in " << profile_sec << " sec; mem = " << worker_mem 
					<< "%" << endl;
			}

			break;

		case STATUS_DATABASE_SUCCESS:

			// Receive the:
			//	- list of accessions
			//	- database index
			//	- time to build database
			//	- %worker memory used
			ptr = mpi_unpack(ptr, acc_list);
			ptr = mpi_unpack(ptr, index);
			ptr = mpi_unpack(ptr, profile_sec);
			ptr = mpi_unpack(ptr, worker_mem);

			for(vector<SraAccession>::const_iterator a = acc_list.begin();a != acc_list.end();++a){

				// Update the accession status in RAM
				m_status[get_accession_index(*a, m_accession_loc)] = m_mpi_status.MPI_TAG;
			}

			// Return the worker to the pool
			m_workers.push_back(m_mpi_status.MPI_SOURCE);

			if(m_verbose){
				cerr << "Worker " << m_mpi_status.MPI_SOURCE << " successfully created database " 
					<< index << " with " << acc_list.size() << " accessions in " << profile_sec << " sec; mem = " 
					<< worker_mem << "%" << endl;
			}

			// Force an immediate write of the status array to disk
			ret = true;

			break;

		case STATUS_DATABASE_FAIL:

			// Receive the:
			//	- list of accessions
			//	- database index
			//	- time to build database
			//	- %worker memory used
			ptr = mpi_unpack(ptr, acc_list);
			ptr = mpi_unpack(ptr, index);
			ptr = mpi_unpack(ptr, profile_sec);
			ptr = mpi_unpack(ptr, worker_mem);

			for(vector<SraAccession>::const_iterator a = acc_list.begin();a != acc_list.end();++a){

				// Update the accession status in RAM
				m_status[get_accession_index(*a, m_accession_loc)] = m_mpi_status.MPI_TAG;
			}

			// Return the worker to the pool
			m_workers.push_back(m_mpi_status.MPI_SOURCE);

			if(m_verbose){
				cerr << "Worker " << m_mpi_status.MPI_SOURCE << " failed to create database " 
					<< index << " with " << acc_list.size() << " accessions in " << profile_sec << " sec; mem = " 
					<< worker_mem << "%" << endl;
			}

			break;

		case STATUS_DATABASE_UPLOAD_FAIL:

			// Receive the:
			//	- list of accessions
			//	- database index
			//	- time to build database
			//	- %worker memory used
			ptr = mpi_unpack(ptr, acc_list);
			ptr = mpi_unpack(ptr, index);
			ptr = mpi_unpack(ptr, profile_sec);
			ptr = mpi_unpack(ptr, worker_mem);

			for(vector<SraAccession>::const_iterator a = acc_list.begin();a != acc_list.end();++a){

				// Update the accession status in RAM
				m_status[get_accession_index(*a, m_accession_loc)] = m_mpi_status.MPI_TAG;
			}

			// Return the worker to the pool
			m_workers.push_back(m_mpi_status.MPI_SOURCE);

			if(m_verbose){
				cerr << "Worker " << m_mpi_status.MPI_SOURCE << " failed to upload database " 
					<< index << " with " << acc_list.size() << " accessions; mem = " 
					<< worker_mem << "%" << endl;
			}

			break;

		default:
			throw __FILE__ ":process_event: Invalid message tag in main event loop";
	};

	delete [] buffer;

	return ret;
}

bool schedule_database_creation(unordered_map< BloomParam, deque<SraAccession> > &m_bloom, 
	size_t &m_database_index, deque<int> &m_workers, 
	const MAP<size_t /*log_2 Bloom length*/, size_t /*number of filters/database file*/> &m_database_size,
	bool m_verbose)
{
	for(unordered_map< BloomParam, deque<SraAccession> >::iterator i = m_bloom.begin();i != m_bloom.end();++i){

		MAP<size_t /*log_2 Bloom length*/, size_t /*number of filters/database file*/>::const_iterator size_iter =
			m_database_size.find(i->first.log_2_filter_len);

		if( size_iter == m_database_size.end() ){
			throw __FILE__ ":schedule_database_creation: Unable to lookup number of filters per database file";
		}

		if(i->second.size() >= size_iter->second){

			if( m_workers.empty() ){

				// We would like to create the database, but we need to wait for a worker to become available
				return true;
			}

			// If the number of filters is zero, this is a special signal to pack all remaining 
			// accessions into a database file. This is needed for the last steps of the SRA database
			// construction.
			const size_t len = (size_iter->second == 0) ? i->second.size() : size_iter->second;

			deque<SraAccession> acc;

			for(size_t j = 0;j < len;++j){

				// Remove an accession from the list of pending accessions 
				// and add it to the active list
				acc.push_back( i->second.front() );
				i->second.pop_front();
			}

			unsigned int buffer_size = 
				mpi_size(m_database_index) + 
				mpi_size(i->first) + 
				mpi_size(acc);

			unsigned char *buffer = new unsigned char[buffer_size];

			if(buffer == NULL){
				throw __FILE__ ":schedule_database_creation: Unable to allocate buffer";
			}

			unsigned char *ptr = buffer;

			ptr = mpi_pack(ptr, m_database_index);
			ptr = mpi_pack(ptr, i->first);
			ptr = mpi_pack(ptr, acc);

			const int w = m_workers.back();
			m_workers.pop_back();

			if( MPI_Send( buffer, buffer_size, MPI_BYTE, w, SCHEDULE_DATABASE, MPI_COMM_WORLD ) != MPI_SUCCESS ){
				throw __FILE__ ":schedule_database_creation: Error sending SCHEDULE_DATABASE to worker";
			}

			delete [] buffer;

			if(m_verbose){
				cerr << "Requested worker " << w << " create database file " << m_database_index 
					<< " with " << acc.size()  << " Bloom filters" << endl;
			}

			// Increment the database file index (to make sure that all database files have unique names)
			++m_database_index;

			// Have we exhuasted all of the accession for this combination of Bloom filter parameters?
			if( i->second.empty() ){

				// Remove elements that do not have at least one associated SRA accession
				m_bloom.erase(i);
			}

			return true;
		}
	}

	return false;
}

bool schedule_bloom_filter_creation(deque<SraAccession> &m_download, 
	deque<int> &m_workers, const vector< pair<SraAccession, size_t /*location*/> > &m_accession_loc, 
	const string &m_metadata_file, bool m_verbose)
{
	// Do we have have already-downloaded SRA records to turn into Bloom filters? Or available workers?
	if( m_download.empty() || m_workers.empty() ){
		return false;
	}

	const SraAccession acc = m_download.front();

	m_download.pop_front();

	// Read the SRA metadata
	FilterInfo info;

	ifstream fin(m_metadata_file.c_str(), ios::binary);

	if(!fin){
		throw __FILE__ ":schedule_bloom_filter_creation: Unable to open metadata file for reading";
	}

	vector< pair<SraAccession, size_t /*location*/> >::const_iterator loc_iter = 
		lower_bound( m_accession_loc.begin(), m_accession_loc.end(), acc, search_by_accession() );

	if( loc_iter == m_accession_loc.end() ){
		throw __FILE__ ":schedule_bloom_filter_creation: Unable to lookup accesion metadata location";
	}
	
	fin.seekg(loc_iter->second, ios_base::beg);

	if(!fin){
		throw __FILE__ ":schedule_bloom_filter_creation: Error seeking position in metadata file for reading";
	}

	binary_read(fin, info);

	// Pack the message
	unsigned int buffer_size = mpi_size(acc) + mpi_size(info);

	unsigned char* buffer = new unsigned char [buffer_size];

	if(buffer == NULL){
		throw __FILE__ ":schedule_bloom_filter_creation: Unable to allocate buffer";
	}

	unsigned char* ptr = buffer;

	// Pack the buffer
	ptr = mpi_pack(ptr, acc);
	ptr = mpi_pack(ptr, info);

	const int w = m_workers.front();

	m_workers.pop_front();

	// Send the buffer to the workers	
	if( MPI_Send( buffer, buffer_size, MPI_BYTE, w, SCHEDULE_BLOOM, MPI_COMM_WORLD ) != MPI_SUCCESS ){
		throw __FILE__ ":schedule_database_creation: Error sending SCHEDULE_BLOOM to worker";
	}

	delete [] buffer;

	if(m_verbose){

		cerr << "Scheduled Bloom filter creation for " << accession_to_str(acc) 
			<< " using worker " << w << endl;
	}

	return true;
}

bool schedule_streaming_bloom_filter_creation(size_t &m_curr_sra, deque<SraAccession> &m_retry, 
	deque<int> &m_workers, unsigned char* m_status, 
	const vector< pair<SraAccession, size_t /*location*/> > &m_accession_loc, 
	const string &m_metadata_file, const unsigned int &m_num_retry, bool m_verbose)
{
	const size_t num_sra = m_accession_loc.size();

	if( m_workers.empty() ){

		// We would like to stream a new SRA record, but can't since there are no free workers.
		return false;
	}

	SraAccession acc = INVALID_ACCESSION;

	// Retries take precedence over unprocessed records
	if( !m_retry.empty() ){

		acc = m_retry.front();

		m_retry.pop_front();
	}
	else{

		while(m_curr_sra < num_sra){

			// Does the current SRA record need to be downloaded and converted into a Bloom filter?
			bool valid_sra = false;

			switch(m_status[m_curr_sra]){
				case STATUS_DOWNLOAD_FAIL:
					break;
				case STATUS_INIT:
				case STATUS_DOWNLOAD_SUCCESS: // When restoring from a non-streaming run, all downloads
				case STATUS_DOWNLOAD_FAIL_1:  // and attempted downloads, will be streamed again.
				case STATUS_DOWNLOAD_FAIL_2:
				case STATUS_DOWNLOAD_FAIL_3:
				case STATUS_DOWNLOAD_FAIL_4:
				case STATUS_DOWNLOAD_FAIL_5:
				case STATUS_DOWNLOAD_FAIL_6:
				case STATUS_DOWNLOAD_FAIL_7:
				case STATUS_DOWNLOAD_FAIL_8:
				case STATUS_DOWNLOAD_FAIL_9:
				case STATUS_DOWNLOAD_FAIL_10:

					valid_sra = true;

					// Mark the status a "failure" in case we never get a response from
					// the worker assigned to download
					m_status[m_curr_sra] = STATUS_BLOOM_FAIL_1;

					break;
				case STATUS_BLOOM_FAIL_1:

					if(m_num_retry > 1){

						valid_sra = true;

						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_BLOOM_FAIL_2;
					}

					break;
				case STATUS_BLOOM_FAIL_2:

					if(m_num_retry > 2){

						valid_sra = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_BLOOM_FAIL_3;
					}

					break;
				case STATUS_BLOOM_FAIL_3:

					if(m_num_retry > 3){

						valid_sra = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_BLOOM_FAIL_4;
					}

					break;
				case STATUS_BLOOM_FAIL_4:

					if(m_num_retry > 4){

						valid_sra = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_BLOOM_FAIL_5;
					}

					break;
				case STATUS_BLOOM_FAIL_5:

					if(m_num_retry > 5){

						valid_sra = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_BLOOM_FAIL_6;
					}

					break;
				case STATUS_BLOOM_FAIL_6:

					if(m_num_retry > 6){

						valid_sra = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_BLOOM_FAIL_7;
					}

					break;
				case STATUS_BLOOM_FAIL_7:

					if(m_num_retry > 7){

						valid_sra = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_BLOOM_FAIL_8;
					}

					break;
				case STATUS_BLOOM_FAIL_8:

					if(m_num_retry > 8){

						valid_sra = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_BLOOM_FAIL_9;
					}

					break;
				case STATUS_BLOOM_FAIL_9:

					if(m_num_retry > 9){

						valid_sra = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_BLOOM_FAIL_10;
					}

					break;
				case STATUS_BLOOM_FAIL_10:

					if(m_num_retry > 10){

						valid_sra = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_BLOOM_FAIL;
					}

					break;
				
				// Enumerate all of the allowed STATUS flags to make sure we throw an error
				// if we do not explicitly handle *every* flag.
				case STATUS_BLOOM_SUCCESS:
				case STATUS_BLOOM_FAIL:
				case STATUS_BLOOM_INVALID:
				case STATUS_DATABASE_SUCCESS:
				case STATUS_DATABASE_FAIL:
				case STATUS_DATABASE_UPLOAD_FAIL:
				case STATUS_SKIPPED:
					break;
				default:
					throw __FILE__ ":schedule_streaming_bloom_filter_creation: Unknown SRA status entry";
			};

			if(valid_sra){
				break;
			}

			++m_curr_sra;
		}

		if(m_curr_sra >= num_sra){

			// There are no more SRA records to stream
			return false;
		}

		acc = m_accession_loc[m_curr_sra].first;
	}

	if(acc == INVALID_ACCESSION){
		throw __FILE__ ":schedule_streaming_bloom_filter_creation: Invalid SRA accession";
	}

	// Read the SRA metadata
	FilterInfo info;

	ifstream fin(m_metadata_file.c_str(), ios::binary);

	if(!fin){
		throw __FILE__ ":schedule_streaming_bloom_filter_creation: Unable to open metadata file for reading";
	}

	vector< pair<SraAccession, size_t /*location*/> >::const_iterator loc_iter = 
		lower_bound( m_accession_loc.begin(), m_accession_loc.end(), acc, search_by_accession() );

	if( loc_iter == m_accession_loc.end() ){
		throw __FILE__ ":schedule_streaming_bloom_filter_creation: Unable to lookup accesion metadata location";
	}
	
	fin.seekg(loc_iter->second, ios_base::beg);

	if(!fin){
		throw __FILE__ ":schedule_streaming_bloom_filter_creation: Error seeking position in metadata file for reading";
	}

	binary_read(fin, info);

	// Pack the message
	unsigned int buffer_size = mpi_size(acc) + mpi_size(info);

	unsigned char* buffer = new unsigned char [buffer_size];

	if(buffer == NULL){
		throw __FILE__ ":schedule_streaming_bloom_filter_creation: Unable to allocate buffer";
	}

	unsigned char* ptr = buffer;

	// Pack the buffer
	ptr = mpi_pack(ptr, acc);
	ptr = mpi_pack(ptr, info);

	const int w = m_workers.front();

	m_workers.pop_front();

	// Send the buffer to the workers	
	if( MPI_Send( buffer, buffer_size, MPI_BYTE, w, SCHEDULE_BLOOM, MPI_COMM_WORLD ) != MPI_SUCCESS ){
		throw __FILE__ ":schedule_streaming_bloom_filter_creation: Error sending SCHEDULE_BLOOM to worker";
	}

	delete [] buffer;

	if(m_verbose){

		cerr << "Scheduled streaming Bloom filter creation for " << accession_to_str(acc) 
			<< " using worker " << w << endl;
	}

	++m_curr_sra;

	return true;
}

bool schedule_download(size_t &m_curr_sra, deque<SraAccession> &m_retry, deque<int> &m_workers, 
	unsigned char* m_status, 
	const vector< pair<SraAccession, size_t> > &m_accession_loc, const unsigned int &m_num_retry, 
	bool m_verbose)
{
	const size_t num_sra = m_accession_loc.size();

	if( m_workers.empty() ){

		// We would like to download a new SRA records, but can't since there are no free workers.
		return false;
	}

	SraAccession acc = INVALID_ACCESSION;

	// Retries take precedence over unprocessed records
	if( !m_retry.empty() ){

		acc = m_retry.front();

		m_retry.pop_front();
	}
	else{

		while(m_curr_sra < num_sra){

			// Does the current SRA record need to be downloaded?
			bool valid_download = false;

			switch(m_status[m_curr_sra]){
				case STATUS_INIT:
				case STATUS_BLOOM_FAIL:
				case STATUS_BLOOM_FAIL_1:	// If we are restoring from a streaming run,
				case STATUS_BLOOM_FAIL_2:	// we need to retry failed Bloom filters
				case STATUS_BLOOM_FAIL_3:
				case STATUS_BLOOM_FAIL_4:
				case STATUS_BLOOM_FAIL_5:
				case STATUS_BLOOM_FAIL_6:
				case STATUS_BLOOM_FAIL_7:
				case STATUS_BLOOM_FAIL_8:
				case STATUS_BLOOM_FAIL_9:
				case STATUS_BLOOM_FAIL_10:

					valid_download = true;

					// Mark the status a "failure" in case we never get a response from
					// the worker assigned to download
					m_status[m_curr_sra] = STATUS_DOWNLOAD_FAIL_1;
					break;
				case STATUS_DOWNLOAD_SUCCESS:
				case STATUS_DOWNLOAD_FAIL:
					break;
				case STATUS_DOWNLOAD_FAIL_1:

					if(m_num_retry > 1){

						valid_download = true;

						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_DOWNLOAD_FAIL_2;
					}

					break;
				case STATUS_DOWNLOAD_FAIL_2:

					if(m_num_retry > 2){

						valid_download = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_DOWNLOAD_FAIL_3;
					}

					break;
				case STATUS_DOWNLOAD_FAIL_3:

					if(m_num_retry > 3){

						valid_download = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_DOWNLOAD_FAIL_4;
					}

					break;
				case STATUS_DOWNLOAD_FAIL_4:

					if(m_num_retry > 4){

						valid_download = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_DOWNLOAD_FAIL_5;
					}

					break;
				case STATUS_DOWNLOAD_FAIL_5:

					if(m_num_retry > 5){

						valid_download = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_DOWNLOAD_FAIL_6;
					}

					break;
				case STATUS_DOWNLOAD_FAIL_6:

					if(m_num_retry > 6){

						valid_download = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_DOWNLOAD_FAIL_7;
					}

					break;
				case STATUS_DOWNLOAD_FAIL_7:

					if(m_num_retry > 7){

						valid_download = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_DOWNLOAD_FAIL_8;
					}

					break;
				case STATUS_DOWNLOAD_FAIL_8:

					if(m_num_retry > 8){

						valid_download = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_DOWNLOAD_FAIL_9;
					}

					break;
				case STATUS_DOWNLOAD_FAIL_9:

					if(m_num_retry > 9){

						valid_download = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_DOWNLOAD_FAIL_10;
					}

					break;
				case STATUS_DOWNLOAD_FAIL_10:

					if(m_num_retry > 10){

						valid_download = true;
						
						// Increment the "failure" status in case we never get a response from
						// the worker assigned to download
						m_status[m_curr_sra] = STATUS_DOWNLOAD_FAIL;
					}

					break;
				
				// Enumerate all of the allowed STATUS flags to make sure we throw an error
				// if we do not explicitly handle *every* flag.
				case STATUS_BLOOM_SUCCESS:
				case STATUS_BLOOM_INVALID:
				case STATUS_DATABASE_SUCCESS:
				case STATUS_DATABASE_FAIL:
				case STATUS_DATABASE_UPLOAD_FAIL:
				case STATUS_SKIPPED:
					break;
				default:
					throw __FILE__ ":schedule_download: Unknown SRA status entry";
			};

			if(valid_download){
				break;
			}

			++m_curr_sra;
		}

		if(m_curr_sra >= num_sra){

			// There are no more SRA records to download
			return false;
		}

		acc = m_accession_loc[m_curr_sra].first;
	}

	if(acc == INVALID_ACCESSION){
		throw __FILE__ ":schedule_download: Invalid SRA accession";
	}

	const int w = m_workers.front();

	m_workers.pop_front();

	unsigned int buffer_size = mpi_size(acc);

	unsigned char* buffer = new unsigned char [buffer_size];

	if(buffer == NULL){
		throw __FILE__ ":schedule_download: Unable to allocate buffer";
	}

	mpi_pack(buffer, acc);

	if( MPI_Send( buffer, buffer_size, MPI_BYTE, w, SCHEDULE_DOWNLOAD, MPI_COMM_WORLD ) != MPI_SUCCESS ){
		throw __FILE__ ":schedule_download: Error sending SCHEDULE_DOWNLOAD to worker";
	}

	delete [] buffer;

	if(m_verbose){

		cerr << "Scheduled worker " << w << " to download " 
			<< accession_to_str(m_accession_loc[m_curr_sra].first) << " ("
			<< (100.0*m_curr_sra)/num_sra << "% complete)" << endl;
	}

	++m_curr_sra;

	return true;
}

vector<string> get_rank_names()
{
	vector<string> ret(mpi_numtasks);

	int name_len = 0;
	char name[MPI_MAX_PROCESSOR_NAME];

	if(MPI_Get_processor_name(name, &name_len) != MPI_SUCCESS){
		throw __FILE__ ":get_rank_names: Unable to obtain processor name";
	}

	ret[0] = name;

	for(int w = 1;w < mpi_numtasks;++w){

		if(MPI_Recv(name, MPI_MAX_PROCESSOR_NAME, MPI_BYTE, w, 
			PROCESSOR_NAME, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS){

			throw __FILE__ ":get_rank_names: Unable to receive name";
		}

		ret[w] = name;
	}

	return ret;
}