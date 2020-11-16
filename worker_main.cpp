#include "maestro.h"
#include "sra_accession.h"
#include "bloom.h"
#include "file_util.h"
#include "mem_usage.h"

#include <sstream>

#include <unistd.h> // unlink
#include <math.h> // ceil

using namespace std;

extern int mpi_rank;
extern int mpi_numtask;

const string prefetch_cmd = "/home/ubuntu/SRA/bin/prefetch";

void send_processor_name();
void process_database_request(unsigned char* m_buffer, const string &m_database_dir, 
	const string &m_bloom_dir, const MaestroOptions &m_opt);
void process_bloom_request(unsigned char* m_buffer, const string &m_bloom_dir, 
	const string &m_sra_dir, const MaestroOptions &m_opt);
void process_download_request(unsigned char* m_buffer, const string &m_sra_dir, 
	const MaestroOptions &m_opt);

void worker_main(MaestroOptions &m_opt)
{
	bool init = true;

	// Send the processor name to the maestro
	send_processor_name();

	// Was rank 0 successfull in initialization?
	broadcast(init, mpi_rank, 0);

	if(!init){
		return;
	}

	string local_sra_repository_dir;

	broadcast(local_sra_repository_dir, mpi_rank, 0);

	// The main event loop
	MPI_Status mpi_status;

	while(true){

		// Check for incomming messages from rank 0
		if(MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status) != MPI_SUCCESS){
			throw __FILE__ ":worker_main: Error in MPI_Probe";
		}

		#ifdef MEM_DEBUG
		cerr << "[" << mpi_rank << "] starting with " << memory_usage() << "% mem" << endl;
		#endif // MEM_DEBUG

		int buffer_size = 0;

		if(MPI_Get_count(&mpi_status, MPI_BYTE, &buffer_size) != MPI_SUCCESS){
			throw __FILE__ ":worker_main: Error obtaining main even buffer size with MPI_Get_count";
		}

		unsigned char *buffer = NULL;
		
		if(buffer_size > 0){
			
			buffer = new unsigned char [buffer_size];

			if(buffer == NULL){
				throw __FILE__ ":worker_main: Unable to allocate main event buffer";
			}
		}

		if(MPI_Recv(buffer, buffer_size, MPI_BYTE, mpi_status.MPI_SOURCE, 
			mpi_status.MPI_TAG, MPI_COMM_WORLD, &mpi_status) != MPI_SUCCESS){

			throw __FILE__ ":worker_main: Unable to receive main even buffer";
		}

		if(mpi_status.MPI_TAG == MAESTRO_QUIT){

			// All done! Since MAESTRO_QUIT is a zero byte message, we do not
			// need to delete the main even message buffer. 
			break;
		}

		switch(mpi_status.MPI_TAG){

			case SCHEDULE_DATABASE:
				process_database_request(buffer, m_opt.scratch_database_dir, m_opt.scratch_bloom_dir, m_opt);
				break;
			case SCHEDULE_BLOOM:
				process_bloom_request(buffer, m_opt.scratch_bloom_dir, local_sra_repository_dir, m_opt);
				break;
			case SCHEDULE_DOWNLOAD:
				process_download_request(buffer, local_sra_repository_dir, m_opt);
				break;
			default:
				throw __FILE__ ":worker_main: Received an illegal message tag in the main even loop";
		}

		if(buffer != NULL){
			delete [] buffer;
		}

		#ifdef MEM_DEBUG
		cerr << "[" << mpi_rank << "] ended with " << memory_usage() << "% mem" << endl;
		#endif // MEM_DEBUG
	}
}

void process_database_request(unsigned char* m_buffer, const string &m_database_dir, 
	const string &m_bloom_dir, const MaestroOptions &m_opt)
{
	// Unpack the:
	// 	database index (size_t)
	//	BloomParam
	//	vector (or deque) of SraAccession
	size_t database_index = 0;
	BloomParam param;
	vector<SraAccession> acc;

	unsigned char *ptr = m_buffer;

	ptr = mpi_unpack(ptr, database_index);
	ptr = mpi_unpack(ptr, param);
	ptr = mpi_unpack(ptr, acc);

	if( acc.empty() ){
		throw __FILE__ ":process_database_request: No accessions!";
	}

	// Note that while we may be packing *fewer* than MAX_NUM_FILTER_CHUNK Bloom filters into a file,
	// we should never be packing more ...
	if(acc.size() > MAX_NUM_FILTER_CHUNK){
		throw __FILE__ ":process_database_request: Exceeded the maximum number of allowed accessions per database file";
	}

	double profile = MPI_Wtime();

	// The database filename
	stringstream database_filename;

	database_filename << "sra." << database_index << ".db";

	const string local_filename = m_database_dir + PATH_SEPARATOR + database_filename.str();

	// The inventory of Bloom filter files to combine into a single database file
	deque<string> bloom_files;

	for(vector<SraAccession>::const_iterator i = acc.begin();i != acc.end();++i){
		bloom_files.push_back(m_bloom_dir + PATH_SEPARATOR + accession_to_str(*i) + ".bloom");
	}

	int ret = STATUS_DATABASE_SUCCESS;

	try{
		if( !build_db(local_filename, param, bloom_files) ){
			ret = STATUS_DATABASE_FAIL;
		}
		else{

			if(!m_opt.save_bloom){

				// We have successfully created a database file. It is safe to delete the Bloom 
				// filter files
				for(deque<string>::const_iterator i = bloom_files.begin();i != bloom_files.end();++i){

					if( unlink( i->c_str() ) != 0 ){
						cerr << "[" << mpi_rank << "] Unable to remove Bloom filter " << *i << endl;
					}
				}
			}
		}
	}
	catch(...){
		ret = STATUS_DATABASE_FAIL;
	}

	if(ret == STATUS_DATABASE_SUCCESS){

		// Success! Now we need to copy to S3
		stringstream script;

		if(m_opt.save_db){

			// Copy the local database file
			script << "aws s3 cp " << local_filename << " " << m_opt.s3_bucket << PATH_SEPARATOR << database_filename.str()
				<< " > /dev/null 2>&1";
		}
		else{
			
			// Move the local database file
			script << "aws s3 mv " << local_filename << " " << m_opt.s3_bucket << PATH_SEPARATOR << database_filename.str()
				<< " > /dev/null 2>&1";
		}

		// If the s3_no_write is true, then we will skip the actual writing of database files to s3.
		// This is useful when debugging.
		if(m_opt.s3_no_write == false){ // <-- double negative!

			if(system( script.str().c_str() ) != 0){
				ret = STATUS_DATABASE_UPLOAD_FAIL;
			}
		}
		else{
		
			// Since we did not move the database file to s3, we need to check to see if the
			// local copy should be saved.
			if(m_opt.save_db == false){
				if( unlink( local_filename.c_str() ) != 0 ){
					cerr << "[" << mpi_rank << "] Unable to remove database file " << local_filename << endl;
				}
			}
		}
	}

	profile = MPI_Wtime() - profile;

	const float worker_mem = memory_usage();

	unsigned int buffer_size = mpi_size(acc) + mpi_size(database_index) + mpi_size(profile) + mpi_size(worker_mem);
	unsigned char* buffer = new unsigned char[buffer_size];

	if(buffer == NULL){
		throw __FILE__ ":process_database_request: Unable to allocate buffer";
	}

	ptr = buffer;

	ptr = mpi_pack(ptr, acc);
	ptr = mpi_pack(ptr, database_index);
	ptr = mpi_pack(ptr, profile);
	ptr = mpi_pack(ptr, worker_mem);

	if( MPI_Send( buffer, buffer_size, MPI_BYTE, 0, ret, MPI_COMM_WORLD ) != MPI_SUCCESS ){
		throw __FILE__ ":process_database_request: Error sending to maestro";
	}

	delete [] buffer;
}

void process_bloom_request(unsigned char* m_buffer, const string &m_bloom_dir, 
	const string &m_sra_dir, const MaestroOptions &m_opt)
{
	unsigned char *ptr = m_buffer;

	SraAccession acc;
	FilterInfo info;

	ptr = mpi_unpack(ptr, acc);
	ptr = mpi_unpack(ptr, info);

	// Make sure that the acc and the info.run_accession are the same
	if(acc != info.run_accession){
		throw __FILE__ ":process_bloom_request: Inconsistent accession and FilterInfo::accession values";
	}

	const string acc_str = accession_to_str(acc);

	double profile = MPI_Wtime();

	unsigned char ret = STATUS_BLOOM_SUCCESS;

	//#define	DOWNLOAD_SRA_BLOOM_STREAM

	#ifdef DOWNLOAD_SRA_BLOOM_STREAM

	// There can be a problem with "prefetch" resuming an interrupted SRA download. To ensure a reliable download, cleanup any
	// partially download SRA files for this accession.
	remove_sra_files(m_sra_dir, acc_str);

	stringstream script;

	#ifdef LOG_PREFETCH
	script << prefetch_cmd << " --max-size " << m_opt.max_sra_file_size_GB << "G " << acc_str 
		<< " >> /scratch/bigsi++/log/prefetch." << mpi_rank << ".txt 2>&1";
	#else

	// Pipe all prefetch output to /dev/null to avoid cluttering up the log files
	script << prefetch_cmd << " --max-size " << m_opt.max_sra_file_size_GB << "G " << acc_str << " > /dev/null 2>&1";
	#endif // LOG_PREFETCH

	if(system( script.str().c_str() ) != 0){
		ret = STATUS_BLOOM_FAIL;
	}

	#endif // DOWNLOAD_SRA_BLOOM_STREAM

	BloomParam param;
	BloomProgress progress;
	
	// If we are download SRA files prior to Bloom filter creation, only proceed if the
	// prefetch command was successful
	if(ret == STATUS_BLOOM_SUCCESS){

		ret = make_bloom_filter(acc, info, param, progress, m_bloom_dir, m_opt);

		// Handle aligned colorspace reads as a special case. According to https://github.com/ncbi/ncbi-vdb/issues/31
		// aligned colorspace SRA records are "broken" and will fail when trying to read primary alignments 
		// followed by unaligned reads. The signature of this failure is that all primary alignments will be
		// read successfully and no unaligned reads will be successfully read.
		if( (ret == STATUS_BLOOM_FAIL) && (progress.num_primary_align > 0) && 
			(progress.curr_primary_align == progress.num_primary_align) &&
			(progress.num_unaligned_read > 0) && (progress.curr_unaligned_read == 0) ){
			
			ret = make_bloom_filter(acc, info, param, progress, m_bloom_dir, m_opt, true /*force unaligned*/);
		}
	}

	if( (ret != STATUS_BLOOM_SUCCESS) || !m_opt.save_sra ){

		// Clean up any partially downloaded "".sra" files along with SRA record-specific helper files (like ".sra.cache" files)
		// Helper files (i.e. .sra.cache) can be downloaded even when *streaming* SRA data.
		remove_sra_files(m_sra_dir, acc_str);
	}

	profile = MPI_Wtime() - profile;

	const float worker_mem = memory_usage();

	unsigned int buffer_size = 0;
	unsigned char* buffer = NULL;

	if(ret == STATUS_BLOOM_SUCCESS){

		// Send the:
		//	- accession
		//	- Bloom param
		//	- Bloom progress
		//  - time to build bloom filter
		// 	- %worker memory used

		buffer_size = mpi_size(acc) + mpi_size(param) + mpi_size(progress)
			+ mpi_size(profile) + mpi_size(worker_mem);

		buffer = new unsigned char[buffer_size];

		if(buffer == NULL){
			throw __FILE__ ":process_bloom_request: Unable to allocate buffer (success)";
		}

		ptr = buffer;

		ptr = mpi_pack(ptr, acc);
		ptr = mpi_pack(ptr, param);
		ptr = mpi_pack(ptr, progress);
		ptr = mpi_pack(ptr, profile);
		ptr = mpi_pack(ptr, worker_mem);
	}
	else{
		
		// Send the:
		//	- accession
		//	- Bloom progress
		//  - time to build bloom filter
		// 	- %worker memory used

		buffer_size = mpi_size(acc) + mpi_size(progress) + mpi_size(profile) + mpi_size(worker_mem);

		buffer = new unsigned char[buffer_size];

		if(buffer == NULL){
			throw __FILE__ ":process_bloom_request: Unable to allocate buffer (failure)";
		}

		ptr = buffer;

		ptr = mpi_pack(ptr, acc);
		ptr = mpi_pack(ptr, progress);
		ptr = mpi_pack(ptr, profile);
		ptr = mpi_pack(ptr, worker_mem);
	}

	if( MPI_Send( buffer, buffer_size, MPI_BYTE, 0, ret, MPI_COMM_WORLD ) != MPI_SUCCESS ){
		throw __FILE__ ":process_bloom_request: Error sending result to maestro";
	}

	delete [] buffer;
}

void process_download_request(unsigned char* m_buffer, const string &m_sra_dir, const MaestroOptions &m_opt)
{
	// Unpack the accession to download from the buffer
	SraAccession acc;
	double profile = MPI_Wtime();

	mpi_unpack(m_buffer, acc);

	const string acc_str = accession_to_str(acc);

	// There can be a problem with "prefetch" resuming an interrupted SRA download. To ensure a reliable download, cleanup any
	// partially download SRA files for this accession.
	remove_sra_files(m_sra_dir, acc_str);

	stringstream script;

	#ifdef LOG_PREFETCH
	script << prefetch_cmd << " --max-size " << m_opt.max_sra_file_size_GB << "G " << acc_str 
		<< " >> /scratch/bigsi++/log/prefetch." << mpi_rank << ".txt 2>&1";
	#else

	// Pipe all prefetch output to /dev/null to avoid cluttering up the log files
	script << prefetch_cmd << " --max-size " << m_opt.max_sra_file_size_GB << "G " << acc_str << " > /dev/null 2>&1";
	#endif // LOG_PREFETCH

	int ret = system( script.str().c_str() );

	// Make sure that we obtained an SRA file (there is at least one example of prefetch downloading a .sra.cache file
	// for DRR000741, but no .sra file)
	if( (ret == 0) && !is_file(m_sra_dir + PATH_SEPARATOR + acc_str + ".sra") ){

		// There is no ".sra" file, so report this download as a failure
		ret = PREFETCH_NO_SRA;
	}

	profile = MPI_Wtime() - profile;

	const float worker_mem = memory_usage();

	unsigned int buffer_size = 0;
	unsigned char *buffer = NULL;

	if(ret == 0){

		buffer_size = mpi_size(acc) + mpi_size(profile) + mpi_size(worker_mem);
		buffer = new unsigned char [buffer_size];

		if(buffer == NULL){
			throw __FILE__ ":process_download_request: Unable to allocate buffer (success)";
		}

		unsigned char* ptr = buffer;

		ptr = mpi_pack(ptr, acc);
		ptr = mpi_pack(ptr, profile);
		ptr = mpi_pack(ptr, worker_mem);
		
		// Success! Return the SRA accession
		if( MPI_Send( buffer, buffer_size, MPI_BYTE, 0, STATUS_DOWNLOAD_SUCCESS, MPI_COMM_WORLD ) != MPI_SUCCESS ){
			throw __FILE__ ":process_download_request: Error sending STATUS_DOWNLOAD_SUCCESS to maestro";
		}
	}
	else{

		buffer_size = mpi_size(acc) + mpi_size(profile) + mpi_size(ret) + mpi_size(worker_mem);
		buffer = new unsigned char [buffer_size];

		if(buffer == NULL){
			throw __FILE__ ":process_download_request: Unable to allocate buffer (failure)";
		}

		unsigned char* ptr = buffer;

		ptr = mpi_pack(ptr, acc);
		ptr = mpi_pack(ptr, profile);
		ptr = mpi_pack(ptr, ret); // <-- return the prefetch error code to help understand the SRA toolkit!
		ptr = mpi_pack(ptr, worker_mem);
		
		// Failure or too large. 
		if( MPI_Send( buffer, buffer_size, MPI_BYTE, 0, STATUS_DOWNLOAD_FAIL, MPI_COMM_WORLD ) != MPI_SUCCESS ){
			throw __FILE__ ":process_download_request: Error sending STATUS_DOWNLOAD_FAIL large status to maestro";
		}

		// Clean up any partially downloaded files (e.g. an orphaned ".sra.cache" file)
		remove_sra_files(m_sra_dir, acc_str);
	}

	if(buffer != NULL){
		delete [] buffer;
	}
}

// Return false if we encounter an error
bool remove_sra_files(const string &m_sra_dir, const string m_accession)
{	
	bool ret = true;

	// Look for files in the local SRA repository that match the 
	// specified accession
	const string search_key = PATH_SEPARATOR + m_accession + ".";

	const size_t num_attempt = 3;

	for(size_t attempt = 0;attempt < num_attempt;++attempt){

		try{
			FindFiles ff(m_sra_dir);

			while( ff.next() ){

				// Find file whose name contains the search_key
				if(ff.name().find(search_key) != string::npos){

					// Both The unlink() function and the bash command "rm -f"
					// occassionally returns an error code. Is the
					// NCBI SRA API keeping one of these files open?
					if(unlink( ff.name().c_str() ) != 0){
						ret = false;
					}

					//const string cmd = "rm -f " + ff.name();

					//if(system( cmd.c_str() ) != 0){
					//	ret = false;
					//}
				}
			}
		}
		catch(...){
			ret = false;
		}

		if(ret){
			return ret;
		}

		sleep(1);
	}

	return ret;
}

void send_processor_name()
{
	int name_len = 0;
	char name[MPI_MAX_PROCESSOR_NAME];

	if(MPI_Get_processor_name(name, &name_len) != MPI_SUCCESS){
		throw __FILE__ ":send_processor_name: Unable to obtain processor name";
	}

	if( MPI_Send( name, MPI_MAX_PROCESSOR_NAME, MPI_BYTE, 0, PROCESSOR_NAME, MPI_COMM_WORLD ) != MPI_SUCCESS ){
		throw __FILE__ ":send_processor_name: Error sending PROCESSOR_NAME to maestro";
	}
}