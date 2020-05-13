#include "sra_stream.h"

#include <iostream>

#include <unistd.h> // For sleep
#include <time.h> // For nanosleep

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

using namespace std;

// When using the VDB API to read SRA data, this is the maximum number 
// of retry attemps before giving up
#define		MAX_RETRY	3

extern int mpi_rank;
extern int mpi_numtasks;

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
	"SRA Download Create Cursor Error"
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
pair<uint64_t, uint64_t> assign_read_range(const uint64_t &m_first_read, const uint64_t &m_num_read);

SRADownloadStatus stream_sra_db_seq(const VDatabase* db_ptr, 
	void per_read_function(const string &m_seq, const unsigned int &m_read_index, const unsigned int &m_read_subindex, void* m_param[]), 
	void* m_param[], StreamStats* m_stat_ptr);

SRADownloadStatus stream_flat_seq(const VTable* tbl_ptr, 
	void per_read_function(const string &m_seq, const unsigned int &m_read_index, const unsigned int &m_read_subindex, void* m_param[]), 
	void* m_param[], StreamStats* m_stat_ptr);

SRADownloadStatus sra_stream(const string &m_accession, 
	void per_read_function(const string &m_seq, const unsigned int &m_read_index, const unsigned int &m_read_subindex, void* m_param[]), 
	void* m_param[], StreamStats* m_stat_ptr /*= NULL*/)
{
	SRADownloadStatus ret = SRADownloadSuccess;

	const VDBManager* mgr = NULL;

	if( VDBManagerMakeRead(&mgr, NULL) ){
		return SRADownloadNetworkFailure;
	}

	// On linux, VDBManagerPathType() appears to work just fine. However, on OS X,
	// when attempting to download an access controlled SRA file (like SRR278685),
	// VDBManagerPathType() hangs indefinately.
	int path_type = VDBManagerPathType( mgr, "%s", m_accession.c_str() );;
	
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

		path_type = VDBManagerPathType( mgr, "%s", m_accession.c_str() );
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
                    		ret = stream_sra_db_seq(db, per_read_function, m_param, m_stat_ptr);
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
				ret = stream_flat_seq(tbl, per_read_function, m_param, m_stat_ptr);
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

SRADownloadStatus stream_sra_db_seq(const VDatabase* db_ptr, 
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
							const pair<uint64_t, uint64_t> read_range = assign_read_range(read_index, num_read);

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

SRADownloadStatus stream_flat_seq(const VTable* tbl_ptr, 
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
					const pair<uint64_t, uint64_t> read_range = assign_read_range(read_index, num_read);

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

pair<uint64_t, uint64_t> assign_read_range(const uint64_t &m_first_read, const uint64_t &m_num_read)
{
	uint64_t chunk = (m_num_read - m_first_read + 1)/mpi_numtasks;
					
	// Each rank is assigned an non-overlapping slice of the SRA
	// file to read. 
	const uint64_t start_index = m_first_read + chunk*mpi_rank;
	
	if( mpi_rank == (mpi_numtasks - 1) ){

		// Read any remainder reads with the last rank
		chunk += (m_num_read - m_first_read + 1)%mpi_numtasks;
	}
	
	const uint64_t stop_index = start_index + chunk;

	return make_pair(start_index, stop_index);
}