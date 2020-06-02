#include "sra_stream.h"
#include "file_util.h"

#include <iostream>

#include <unistd.h> // For sleep
#include <time.h> // For nanosleep

#ifdef _OPENMP
#include <omp.h> // For omp_get_num_threads()
#endif // _OPENMP

// For the num_gen datastructure and helper functions
#include <klib/num-gen.h>

#include <vdb/manager.h>

// For VDBManagerPathType return values (kptDatabase, kptPrereleaseTbl, kptTable)
#include <kdb/manager.h>

// For VDatabase
#include <vdb/database.h>

// Needed for VCursor
#include <vdb/cursor.h>

// For VDatabaseOpenTableRead()
#include <vdb/table.h>

// For openReadCollection (when reading from a file)
#include <ncbi-vdb/NGS.hpp>

using namespace std;

// When using the VDB API to read SRA data, this is the maximum number 
// of retry attemps before giving up
#define		MAX_RETRY	3

// Must match the order enum SRADownloadStatus
const char* SRADownloadErrorStr[] = {
	"Success",
	"SRA Download Network Failure",
	"SRA Download Controlled Access",
	"SRA Download VDB Error",
	"SRA Download List Table Error",
	"SRA Download Read Length Error",
	"SRA Download Cell Data Error",
	"SRA Download Name List Error",
	"SRA Download Add Column Read Error",
	"SRA Download Add Column Read Len Error",
	"SRA Download Cursor Open Error",
	"SRA Download Read Format Error",
	"SRA Download Create Cursor Error",
	"SRA Download Dir Error",
	"SRA Download File Read Error"
};

#define	KB	1024
#define	MB	(KB*KB)
#define	GB	(KB*MB)

// The default "cursor cache" (buffer size?) used by fasterq-dump.c is 5MB
// Can this value be tuned to improve download speed? <-- Nope! This buffer
// speed does not seem influence download speed (since we are accessing each
// read just once). Making this cache smaller does seem to reduce the memory
// footprint (suggesting that it is caching reads in case the user requests
// the same read more than once). By reducing the cache size, we may be able to
// run on compute nodes with smaller amounts of RAM.
#define DEFAULT_CURSOR_CACHE (0 * MB)

// Helper functions cribbed from cmn_iter.c
static bool contains( VNamelist * tables, const char * table );
pair<uint64_t, uint64_t> assign_read_range(const uint64_t &m_first_read, const uint64_t &m_num_read,
	const int &m_rank, const int &m_numtasks);

SRADownloadStatus stream_sra_db_seq(const VDatabase* db_ptr, const int &m_rank, const int &m_numtasks,
	void per_read_function(const string &m_seq, const unsigned int &m_read_index, const unsigned int &m_read_subindex, void* m_param[]), 
	void* m_param[], StreamStats* m_stat_ptr);

SRADownloadStatus stream_flat_seq(const VTable* tbl_ptr, const int &m_rank, const int &m_numtasks,
	void per_read_function(const string &m_seq, const unsigned int &m_read_index, const unsigned int &m_read_subindex, void* m_param[]), 
	void* m_param[], StreamStats* m_stat_ptr);

SRADownloadStatus stream_file(const string &m_accession, const int &m_rank, const int &m_numtasks,
	void per_read_function(const string &m_seq, const unsigned int &m_read_index, const unsigned int &m_read_subindex, void* m_param[]), 
	void* m_param[], StreamStats* m_stat_ptr);

string path_type_to_str(int m_type);

SRADownloadStatus sra_stream(const string &m_accession, const int &m_rank, const int &m_numtasks,
	void per_read_function(const string &m_seq, const unsigned int &m_read_index, const unsigned int &m_read_subindex, void* m_param[]), 
	void* m_param[], StreamStats* m_stat_ptr /*= NULL*/)
{
	SRADownloadStatus ret = SRADownloadSuccess;

	// If m_accession is a valid filesystem path (i.e. a file or a directory), attempt to read the 
	// SRA data from the filesystem
	if( is_path(m_accession) ){
		return stream_file(m_accession, m_rank, m_numtasks, per_read_function, m_param, m_stat_ptr);
	}

	const VDBManager* mgr = NULL;

	if( VDBManagerMakeRead(&mgr, NULL) ){
		return SRADownloadNetworkFailure;
	}

	// On linux, VDBManagerPathType() appears to work just fine. However, on OS X,
	// when attempting to download an access controlled SRA file (like SRR278685),
	// VDBManagerPathType() hangs indefinately.
	int path_type = VDBManagerPathType( mgr, "%s", m_accession.c_str() ) & ~kptAlias;

	// When many processes are trying to access the same SRA record at the
	// same time, VDBManagerPathType will sometimes spuriously return kptNotFound.
	// If this happens, try waiting a fixed amount of time and retrying
	for(unsigned int retry = 0; (path_type == kptNotFound) && (retry < MAX_RETRY);++retry){

		// By waiting a fixed amount of time, we may be creating a "second wave"
		// of ranks that are all calling VDBManagerPathType at the same time
		// (and will cause additiona timeouts). Need to test to see if a random
		// wait time is needed ...
		struct timespec delay;

		delay.tv_sec = 3;
		delay.tv_nsec = 0;

		nanosleep(&delay, NULL);

		path_type = VDBManagerPathType( mgr, "%s", m_accession.c_str() ) & ~kptAlias;
	}

	const VDatabase* db = NULL;
	KNamelist* k_tables = NULL;
	VNamelist* tables = NULL;
	const VTable* tbl = NULL;

	switch(path_type)
	{
		case kptDatabase: 

    		if( VDBManagerOpenDBRead( mgr, &db, NULL, "%s", m_accession.c_str() ) ){
				ret = SRADownloadNetworkFailure;
			}
			else{
				if( VDatabaseListTbl(db, &k_tables) ){
					ret = SRADownloadListTableError;
				}
				else{
					if( VNamelistFromKNamelist(&tables, k_tables) ){
						ret = SRADownloadNameListError;
					}
					else{
						if ( contains( tables, "SEQUENCE" ) ){

							// Ignore PRIMARY_ALIGNMENT and SECONDARY_ALIGNMENT tables. All we care about
							// is sequence data!
                    		ret = stream_sra_db_seq(db, m_rank, m_numtasks, per_read_function, m_param, m_stat_ptr);
                		}
						else{
							// Did not find a "SEQUENCE" table!
							ret = SRADownloadVDBError;
						}

						#ifdef DUMP_TABLES
						{
							uint32_t number_of_table_names = 0;

							VNameListCount(tables, &number_of_table_names);
							cerr << "number_of_table_names = " << number_of_table_names << endl;
							for(uint32_t i = 0;i < number_of_table_names;++i){

								const char* name = NULL;

								VNameListGet(tables, i, &name);

								cerr << "\tname[" << i << "] = " << name << endl;
							}
						}
						#endif // DUMP_TABLES
					}
				}
			}
			
			KNamelistRelease(k_tables);
			VDatabaseRelease(db);

			break;

		case kptPrereleaseTbl:
		case kptTable:
			
			if( VDBManagerOpenTableRead( mgr, &tbl, NULL, "%s", m_accession.c_str() ) ){
				ret = SRADownloadNetworkFailure;
			}
			else{
				ret = stream_flat_seq(tbl, m_rank, m_numtasks, per_read_function, m_param, m_stat_ptr);
			}

			VTableRelease(tbl);

			break;
		default:

			// Invalid VDB path type. This is most likely an unauthorized SRA record
			ret = SRADownloadControlledAccess;
	};

	VDBManagerRelease(mgr);

	return ret;
}

// From the SRA faster-dump.c code
static bool contains( VNamelist * tables, const char * table )
{
    uint32_t found = 0;
    
    return (VNamelistIndexOf(tables, table, &found) == 0);
}

SRADownloadStatus stream_sra_db_seq(const VDatabase* db_ptr, const int &m_rank, const int &m_numtasks,
	void per_read_function(const string &m_seq, const unsigned int &m_read_index, const unsigned int &m_read_subindex, void* m_param[]), 
	void* m_param[], StreamStats* m_stat_ptr)
{
	SRADownloadStatus ret = SRADownloadSuccess;
	
	const VCursor* cur = NULL;
	const VTable* tbl = NULL;
	
	if( VDatabaseOpenTableRead(db_ptr, &tbl, "%s", "SEQUENCE") ){
		ret = SRADownloadNetworkFailure;
	}
	else{
		if( VTableCreateCachedCursorRead(tbl, &cur, DEFAULT_CURSOR_CACHE) ){
			ret = SRADownloadCreateCursorError;
		}
		else{
			uint32_t read_col_id = 0;
			
			if( VCursorAddColumn(cur, &read_col_id, "READ") ){
				ret = SRADownloadAddColumnReadError;
			}
			else{
				uint32_t read_len_col_id = 0;
				
				if( VCursorAddColumn(cur, &read_len_col_id, "READ_LEN") ){
					ret = SRADownloadAddColumnReadLenError;
				}
				else{
				
					if( VCursorOpen(cur) ){
						ret = SRADownloadCursorOpenError;
					}
					else{

						int64_t read_index = 0;
						uint64_t num_read = 0;

						if( VCursorIdRange(cur, read_col_id, &read_index, &num_read) || (read_index < 0) ){
							ret = SRADownloadVDBError;
						}
						else{

							// Each MPI rank gets is assigned a non-overlapping chunk of reads
							const pair<uint64_t, uint64_t> read_range = 
								assign_read_range(read_index, num_read, m_rank, m_numtasks);

							// DEBUG
							//cerr << "[" << mpi_rank << "] stream_sra_db_seq: [" << read_range.first << ", " << read_range.second << "]" << endl;
							//cerr << "[" << mpi_rank << "] read_index = " << read_index << endl;
							//cerr << "[" << mpi_rank << "] num_read = " << num_read << endl;

							uint32_t elem_bits, boff, num_read_len;
							uint32_t *read_len_array = NULL;
							
							String str;

							// Note that the read_index is 1's based (but have taken this into account in read_range)
							for(uint64_t index = read_range.first;index < read_range.second;++index){

								if( VCursorCellDataDirect(cur, index, read_col_id, &elem_bits,
									(const void **)&str.addr, &boff, &str.len ) ){
									
									// The call to VCursorCellDataDirect did not successfully complete,
									// try a few more times before giving up
									struct timespec delay;

									delay.tv_sec = 1;
									delay.tv_nsec = 0;

									unsigned int retry = 0;

									for(;retry < MAX_RETRY;++retry){

										nanosleep(&delay, NULL);

										if( !VCursorCellDataDirect(cur, index, read_col_id, &elem_bits,
											(const void **)&str.addr, &boff, &str.len) ){
											
											// Success!
											break;
										}
									}

									if(retry >= MAX_RETRY){

										ret = SRADownloadCellDataError;
										break;
									}
								}

								if ( (elem_bits != 8) || (boff != 0) ){
								
									ret = SRADownloadReadFormatError;
									break;
								}

								// Only true when storing 8 bit strings
								str.size = str.len;
																
								if( VCursorCellDataDirect( cur, index, read_len_col_id, &elem_bits,
									(const void **)&read_len_array, &boff, &num_read_len) ){
									
									// The call to VCursorCellDataDirect did not successfully complete,
									// try a few more times before giving up
									struct timespec delay;

									delay.tv_sec = 1;
									delay.tv_nsec = 0;

									unsigned int retry = 0;

									for(;retry < MAX_RETRY;++retry){

										nanosleep(&delay, NULL);

										if( !VCursorCellDataDirect( cur, index, read_len_col_id, &elem_bits,
											(const void **)&read_len_array, &boff, &num_read_len) ){
											
											// Success!
											break;
										}
									}

									if(retry >= MAX_RETRY){

										ret = SRADownloadCellDataError;
										break;
									}
								}

								// SRR7841648 is an example of a PacBio dataset with num_read_len == 3
								// This same record (SRR7841648) also has a read with 0 bp ...
								if ( (elem_bits != 32) || (boff != 0) ){

									ret = SRADownloadReadFormatError;
									break;
								}
								
								const char* seq_ptr = str.addr;
		
								for(uint32_t i = 0;i < num_read_len;++i){
									
									per_read_function( string(seq_ptr, read_len_array[i]), 
										index, i + 1 /*sub-read*/, m_param);
									
									if(m_stat_ptr != NULL){

										++m_stat_ptr->num_reads;
										m_stat_ptr->num_bases += read_len_array[i];
									}

									seq_ptr += read_len_array[i];
								}

								// Did we read the expected number of bases?
								//if( (seq_ptr - str.addr) != str.len){
								
									/////////////////////////////////////////////////////////////////////////////////
									// Update -- it appears from the fastq_iter.c file (lines 374 - 391) that we can
									// skip this test! In the fastq_iter.c code, they test for this occurance and,
									// when it happens, they set str.len = sum of read_len_array elements!
									/////////////////////////////////////////////////////////////////////////////////

									//cerr << ":stream_sra_db_seq: Warning! Split read lengths do not equal concatinated read length (sra db)" << endl;
									//cerr << "Expected " << str.len << " bytes, but read " << (seq_ptr - str.addr) << " bytes" << endl;
									//cerr << "Read index = " << index << endl;
									//cerr << "num_read = " << num_read << endl;
									//cerr << "num_read_len = " << num_read_len << endl;

									//for(uint32_t i = 0;i < num_read_len;++i){
									//	cerr << "\tread len[" << i << "] = " << read_len_array[i] << endl;
									//}

									//return SRADownloadReadLengthError;
								//}

								// Don't attempt to deallocate the string returned by VCursorCellDataDirect,
								// doing so generates the error: "pointer being freed was not allocated"
								// StringWhack(&str);
							}
						}
					}
				}
			}
		}
	}
	
	VCursorRelease(cur);
	VTableRelease(tbl);
	
	return ret;
}

SRADownloadStatus stream_flat_seq(const VTable* tbl_ptr, const int &m_rank, const int &m_numtasks,
	void per_read_function(const string &m_seq, const unsigned int &m_read_index, const unsigned int &m_read_subindex, void* m_param[]), 
	void* m_param[], StreamStats* m_stat_ptr)
{
	SRADownloadStatus ret = SRADownloadSuccess;
	
	const VCursor* cur = NULL;

	if( VTableCreateCachedCursorRead(tbl_ptr, &cur, DEFAULT_CURSOR_CACHE) ){
		ret = SRADownloadCreateCursorError;
	}
	else{
	
		uint32_t read_col_id = 0;
		
		if( VCursorAddColumn(cur, &read_col_id, "READ") ){
			ret = SRADownloadAddColumnReadError;
		}
		else{

			if( VCursorOpen(cur) ){
				ret = SRADownloadCursorOpenError;
			}
			else{
				int64_t read_index = 0;
				uint64_t num_read = 0;

				if( VCursorIdRange(cur, read_col_id, &read_index, &num_read) || (read_index < 0) ){
					ret = SRADownloadVDBError;
				}
				else{

					// Each MPI rank gets is assigned a non-overlapping chunk of reads
					const pair<uint64_t, uint64_t> read_range = 
						assign_read_range(read_index, num_read, m_rank, m_numtasks);

					// DEBUG
					//cerr << "[" << mpi_rank << "] stream_flat_seq: [" << read_range.first << ", " << read_range.second << "]" << endl;
					//cerr << "[" << mpi_rank << "] read_index = " << read_index << endl;
					//cerr << "[" << mpi_rank << "] num_read = " << num_read << endl;
					
					uint32_t elem_bits, boff;
					
					String str;

					// Note that the read_index is 1's based (but have taken this into account in read_range)
					for(uint64_t index = read_range.first;index < read_range.second;++index){

						if( VCursorCellDataDirect(cur, index, read_col_id, &elem_bits,
							(const void **)&str.addr, &boff, &str.len ) ){
							
							// The call to VCursorCellDataDirect did not successfully complete,
							// try a few more times before giving up
							struct timespec delay;

							delay.tv_sec = 1;
							delay.tv_nsec = 0;

							unsigned int retry = 0;

							for(;retry < MAX_RETRY;++retry){

								nanosleep(&delay, NULL);

								if( !VCursorCellDataDirect(cur, index, read_col_id, &elem_bits,
									(const void **)&str.addr, &boff, &str.len) ){
									
									// Success!
									break;
								}
							}

							if(retry >= MAX_RETRY){

								ret = SRADownloadCellDataError;
								break;
							}
						}

						if ( (elem_bits != 8) || (boff != 0) ){
							
							ret = SRADownloadReadFormatError;
							break;
						}

						// Only true when storing 8 bit strings
						str.size = str.len;
						
						per_read_function( string(str.addr, str.len), index, 0 /*no sub-read*/, m_param);
						
						if(m_stat_ptr != NULL){

							++m_stat_ptr->num_reads;
							m_stat_ptr->num_bases += str.len;
						}
						
						// Don't attempt to deallocate the string returned by VCursorCellDataDirect,
						// doing so generates the error: "pointer being freed was not allocated"
						// StringWhack(&str);
					}
				}
			}
		}
	}
	
	VCursorRelease(cur);
	
	return ret;
}

pair<uint64_t, uint64_t> assign_read_range(const uint64_t &m_first_read, const uint64_t &m_num_read,
	const int &m_rank, const int &m_numtasks)
{
	uint64_t chunk = (m_num_read - m_first_read + 1)/m_numtasks;
					
	// Each rank is assigned an non-overlapping slice of the SRA
	// file to read. 
	const uint64_t start_index = m_first_read + chunk*m_rank;
	
	if( m_rank == (m_numtasks - 1) ){

		// Read any remainder reads with the last rank
		chunk += (m_num_read - m_first_read + 1)%m_numtasks;
	}
	
	const uint64_t stop_index = start_index + chunk;

	return make_pair(start_index, stop_index);
}

string path_type_to_str(int m_type)
{
	switch(m_type & ~kptAlias){
		case kptNotFound:
			return "kptNotFound";
		case kptBadPath:
			return "kptBadPath";
		case kptFile:
			return "kptFile";
		case kptDir:
			return "kptDir";
		case kptCharDev:
			return "kptCharDev";
		case kptBlockDev:
			return "kptBlockDev";
		case kptFIFO:
			return "kptFIFO";
		case kptZombieFile:
			return "kptZombieFile";
		case kptDataset:
			return "kptDataset";
		case kptDatatype:
			return "kptDatatype";
		case kptDatabase:
			return "kptDatabase";
		case kptTable:
			return "kptTable";
		case kptIndex:
			return "kptIndex";
		case kptColumn:
			return "kptColumn";
		case kptMetadata:
			return "kptMetadata";
		case kptPrereleaseTbl:
			return "kptPrereleaseTbl";
	};

	return "Unknown path type!";
}

SRADownloadStatus stream_file(const string &m_accession, const int &m_rank, const int &m_numtasks,
	void per_read_function(const string &m_seq, const unsigned int &m_read_index, const unsigned int &m_read_subindex, void* m_param[]), 
	void* m_param[], StreamStats* m_stat_ptr)
{
	SRADownloadStatus ret = SRADownloadSuccess;
	
	string sra_dir;
	string sra_file;

	if( is_dir(m_accession) ){

		sra_dir = m_accession;

		// This is the file name that is *local* to sra_dir (since we will chdir to sra_dir)
		sra_file = leaf_path_name(m_accession) + ".sra";
	}
	else{

		sra_dir = parent_dir(m_accession);

		// This is the file name that is local to sra_dir (since we will chdir to sra_dir)
		sra_file = leaf_path_name(m_accession);
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
		return SRADownloadDirError;
	}

	if( chdir( sra_dir.c_str() ) != 0){
		return SRADownloadDirError;
	}

	try{
		// Digest the current sequence group. Trial and error testing on beagle.lanl.gov
		// indicates that 5 reading threads provides the best performance.
		#pragma omp parallel // num_threads(5)
		{
			#ifdef _OPENMP
			const size_t num_thread = omp_get_num_threads();
			const size_t tid = omp_get_thread_num();
			#else
			const size_t num_thread = 1;
			const size_t tid = 0;
			#endif // _OPENMP
			
			StreamStats local_stats;

			ngs::ReadCollection run(  ncbi::NGS::openReadCollection(sra_file) );
			
			// Note that num_read is the number of either paired or
			// unpaired reads. For paired reads, this is half the
			// the number of sequences!
			const size_t num_read = run.getReadCount(ngs::Read::all);
			
			// Each MPI rank gets is assigned a non-overlapping chunk of reads.
			// The returned range is 1's based (to match the SRA API)
			const pair<uint64_t, uint64_t> read_range = 
				assign_read_range(1 /*starting read is 1*/, num_read, m_rank, m_numtasks);

			// For the chunk of reads assigned to this MPI rank, we will *further* divide
			// the chunk amoung all of the threads on the MPI rank.
			const size_t num_local_read = read_range.second - read_range.first;

			size_t chunk = max(size_t(1), num_local_read/num_thread);
										
			// Each thread is assigned an non-overlapping slice of the SRA
			// file to read.
			const size_t start = read_range.first + chunk*tid;
			
			if( tid == (num_thread - 1) ){

				// Assign any "remainder" reads to the last thread
				chunk += num_local_read%num_thread;
			}
							
			ngs::ReadIterator run_iter = 
				ngs::ReadIterator( run.getReadRange ( start, chunk, ngs::Read::all ) );
			
			// Start read counting at zero, since we will be adding this to "start"
			size_t read_count = 0;
						
			while( run_iter.nextRead() ){
				
				// Use 1's based sub-read counting
				size_t sub_read_count = 1;

				while( run_iter.nextFragment() ){

					const string seq = run_iter.getFragmentBases().toString();

					local_stats.num_bases += seq.size();
					++local_stats.num_reads;

					per_read_function(seq, 
						start + read_count, sub_read_count, m_param);

					++sub_read_count;
				}

				++read_count;
			}

			#pragma omp critical
			if(m_stat_ptr != NULL){

				m_stat_ptr->num_reads += local_stats.num_reads;
				m_stat_ptr->num_bases += local_stats.num_bases;
			}	
		}
	}
	catch(...){
		return SRADownloadFileReadError;
	}

	// Return to the original directory after reading each SRA file
	if(chdir(orig_dir) != 0){
		return SRADownloadDirError;
	}
	
	return ret;
}