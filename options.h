#ifndef __OPTIONS
#define	__OPTIONS

#include <unordered_map>
#include <unordered_set>
#include <deque>
#include "mpi_util.h"
#include "hash.h"
#include "date.h"
#include "sra_accession.h"

// We need to protect any commas that appear in template variables
// if we are going to combine them with X Macros
#define SINGLE_ARG(...) __VA_ARGS__

struct SearchOptions
{

	// Allowed output file formats
	enum {
		OUTPUT_CSV,
		OUTPUT_JSON
	};
	
	// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
	// ensure that structure variable are correctly serialized
	#define SEARCH_OPTION_MEMBERS \
		VARIABLE(std::deque<std::string>, query_files) \
		VARIABLE(std::deque<std::string>, query_seq) \
		VARIABLE(std::deque<std::string>, subject_files) \
		VARIABLE(std::string, output_file) \
		VARIABLE(float, threshold) \
		VARIABLE(unsigned int, output_format) \
		VARIABLE(bool, quit)

	#define VARIABLE(A, B) A B;
		SEARCH_OPTION_MEMBERS
	#undef VARIABLE

	SearchOptions(int argc, char* argv[]);
};

struct InventoryOptions
{
	// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
	// ensure that structure variable are correctly serialized
	#define INVENTORY_OPTION_MEMBERS \
		VARIABLE(std::string, metadata_file) \
		VARIABLE(std::string, output_file) \
		VARIABLE(std::unordered_set<std::string>, required_strategy) \
		VARIABLE(std::unordered_set<std::string>, required_source) \
		VARIABLE(std::deque<SraAccession>, include_accessions) \
		VARIABLE(Date, begin_date) \
		VARIABLE(Date, end_date) \
		VARIABLE(bool, list_only) \
		VARIABLE(bool, quit)

	#define VARIABLE(A, B) A B;
		INVENTORY_OPTION_MEMBERS
	#undef VARIABLE

	InventoryOptions()
	{
	};
	
	InventoryOptions(int argc, char* argv[])
	{
		load(argc, argv);
	};
	
	void load(int argc, char* argv[]);
};

struct MaestroOptions
{

	// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
	// ensure that structure variable are correctly serialized
	#define MAESTRO_OPTION_MEMBERS \
		VARIABLE(std::string, metadata_file) \
		VARIABLE(std::string, scratch_bloom_dir) \
		VARIABLE(std::string, scratch_database_dir) \
		VARIABLE(std::string, status_file) \
		VARIABLE(std::string, s3_bucket) \
		VARIABLE(std::deque<SraAccession>, skip_sra) \
		VARIABLE(float, false_positive_probability) \
		VARIABLE(float, download_delay) \
		VARIABLE(unsigned int, min_kmer_count) \
		VARIABLE(unsigned int, kmer_len) \
		VARIABLE(unsigned int, min_log_2_filter_len) \
		VARIABLE(unsigned int, max_log_2_filter_len) \
		VARIABLE(unsigned int, max_sra_file_size_GB) \
		VARIABLE(HashFunction, hash_func) \
		VARIABLE(unsigned int, num_download_attempt) \
		VARIABLE(unsigned int, limit_num_download) \
		VARIABLE(bool, retry_bloom) \
		VARIABLE(bool, save_bloom) \
		VARIABLE(bool, save_db) \
		VARIABLE(bool, save_sra) \
		VARIABLE(bool, s3_no_write) \
		VARIABLE(bool, stream_sra) \
		VARIABLE(bool, verbose) \
		VARIABLE(bool, quit)

	#define VARIABLE(A, B) A B;
		MAESTRO_OPTION_MEMBERS
	#undef VARIABLE

	MaestroOptions()
	{
	};
	
	MaestroOptions(int argc, char* argv[])
	{
		load(argc, argv);
	};
	
	void load(int argc, char* argv[]);
	
	template<class T> friend size_t mpi_size(const T &m_obj);
	template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, 
		T &m_obj);
	template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr,
		const T &m_obj);
};

template<> size_t mpi_size(const MaestroOptions &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const MaestroOptions &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, MaestroOptions &m_obj);

// Human-readable sizes
#define	KB	1024ULL // Don't overflow an int for GB-scale values!
#define	MB	(KB*KB)
#define	GB	(KB*MB)

// The maximum number of groups (i.e. Bloom filters) per file
#define		MAX_NUM_FILTER_CHUNK					2048UL
#define		MAX_DATABASE_FILE_SIZE_IN_GB			64 // In GB
#define		MIN_NUM_FILTERS_FOR_COMPRESSION			64 // Require this many Bloom filters in a slice before compressing
#define		DEFAULT_FALSE_POSITIVE_PROBABILITY		0.25f // Per k-mer
#define		DEFAULT_DOWNLOAD_DELAY					0.0f // Minimum time interval (in seconds) between download/streaming requests
#define		DEFAULT_SUBFILTER_BITS					1024
#define		DEFAULT_COMPRESSION_THRESHOLD			0.75
#define		DEFAULT_KMER_LENGTH						31

#define		MAX_SRA_MIN_KMER_COUNT					15	// This is the MAX_COUNT value from make_bloom.cpp
#define		DEFAULT_SRA_MIN_KMER_COUNT				5
#define		DEFAULT_SEARCH_THRESHOLD 				1.0f
#define		DEFAULT_SEARCH_OUTPUT 					SearchOptions::OUTPUT_JSON
#define		DEFAULT_DOWNLOAD_ATTEMPT				3
#define		DEFAULT_HASH_FUNCTION					MURMUR_HASH32
#define		DEFAULT_MAX_BACKLOG						25
#define		DEFAULT_DOWNLOADS_PER_COMPUTER			1
#define		DEFAULT_MIN_LOG_2_FILTER_LEN			18
#define		DEFAULT_MAX_LOG_2_FILTER_LEN			32 // <-- Must be compatible with the hash function
#define		DEFAULT_MAX_SRA_FILE_SIZE				30 // <-- Max SRA file size in GB
#define		DEFAULT_STATUS_FILE						"./__sra_db_status.bin"

std::ostream& operator<<(std::ostream &m_out, const MaestroOptions &m_opt);

#endif // __OPTIONS
