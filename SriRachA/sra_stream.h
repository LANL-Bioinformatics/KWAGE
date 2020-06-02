#ifndef __SRA_STREAM
#define __SRA_STREAM

#include <string>

typedef enum {
	SRADownloadSuccess, // Must allways have the smallest value to allow max() to find failures
	SRADownloadNetworkFailure,
	SRADownloadControlledAccess,
	SRADownloadVDBError,
	SRADownloadListTableError,
	SRADownloadReadLengthError,
	SRADownloadCellDataError,
	SRADownloadNameListError,
	SRADownloadAddColumnReadError,
	SRADownloadAddColumnReadLenError,
	SRADownloadCursorOpenError,
	SRADownloadReadFormatError,
	SRADownloadCreateCursorError,
	SRADownloadDirError,
	SRADownloadFileReadError
} SRADownloadStatus;

extern const char* SRADownloadErrorStr[];

struct StreamStats
{
	size_t num_reads;
	size_t num_bases;

	StreamStats() :
		num_reads(0), num_bases(0)
	{
		// Do nothing!
	};
};

SRADownloadStatus sra_stream(const std::string &m_accession, const int &m_rank, const int &m_numtasks,
	void per_read_function(const std::string &m_seq, 
		const unsigned int &m_read_index, const unsigned int &m_read_subindex, void* m_param[]), 
	void** m_param, StreamStats* m_stat_ptr = NULL);

#endif // __SRA_STREAM