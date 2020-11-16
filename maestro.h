#ifndef __MAESTO
#define __MAESTO

#include "options.h"
#include "bloom.h"

// The allowed status values and MPI message tags for each SRA (allowed values are from 0 to 255)
#define		STATUS_INIT					1		// Unprocessed
#define		STATUS_DOWNLOAD_SUCCESS		2		// Downloaded

#define		STATUS_DOWNLOAD_FAIL		3		// Unable to download (gave up)
#define		STATUS_DOWNLOAD_FAIL_1		4		// Unable to download (1 attempt)
#define		STATUS_DOWNLOAD_FAIL_2		5		// Unable to download (2 attempts)
#define		STATUS_DOWNLOAD_FAIL_3		6		// Unable to download (3 attempts)
#define		STATUS_DOWNLOAD_FAIL_4		7		// Unable to download (4 attempts)
#define		STATUS_DOWNLOAD_FAIL_5		8		// Unable to download (5 attempts)
#define		STATUS_DOWNLOAD_FAIL_6		9		// Unable to download (6 attempts)
#define		STATUS_DOWNLOAD_FAIL_7		10		// Unable to download (7 attempts)
#define		STATUS_DOWNLOAD_FAIL_8		11		// Unable to download (8 attempts)
#define		STATUS_DOWNLOAD_FAIL_9		12		// Unable to download (9 attempts)
#define		STATUS_DOWNLOAD_FAIL_10		13		// Unable to download (10 attempts)

#define		MAX_NUM_DOWNLOAD_FAIL	STATUS_DOWNLOAD_FAIL_10

#define		STATUS_BLOOM_SUCCESS		14		// Bloom filter file created
#define		STATUS_BLOOM_FAIL			15		// Failed to create a Bloom filter file
#define		STATUS_BLOOM_INVALID		16		// Could not satisfy false positive rate or no data for Bloom filter file
#define		STATUS_DATABASE_SUCCESS		17		// Successfully added to a database file
#define		STATUS_DATABASE_FAIL		18		// Unable to create a database file
#define		STATUS_DATABASE_UPLOAD_FAIL	19		// Failed to upload a database file to S3 storage

// These are used when streaming, which combines downloading and Bloom filter creation.
// This will allow us to retry Bloom filter creation when we encounter a streaming faliure.
#define		STATUS_BLOOM_FAIL_1			20		// Failed to create a Bloom filter file (1 attempt)
#define		STATUS_BLOOM_FAIL_2			21		// Failed to create a Bloom filter file (2 attempts)
#define		STATUS_BLOOM_FAIL_3			22		// Failed to create a Bloom filter file (3 attempts)
#define		STATUS_BLOOM_FAIL_4			23		// Failed to create a Bloom filter file (4 attempts)
#define		STATUS_BLOOM_FAIL_5			24		// Failed to create a Bloom filter file (5 attempts)
#define		STATUS_BLOOM_FAIL_6			25		// Failed to create a Bloom filter file (6 attempts)
#define		STATUS_BLOOM_FAIL_7			26		// Failed to create a Bloom filter file (7 attempts)
#define		STATUS_BLOOM_FAIL_8			27		// Failed to create a Bloom filter file (8 attempts)
#define		STATUS_BLOOM_FAIL_9			28		// Failed to create a Bloom filter file (9 attempts)
#define		STATUS_BLOOM_FAIL_10		29		// Failed to create a Bloom filter file (10 attempts)

#define		STATUS_SKIPPED				30		// Skipped at user request

// Additional MPI message tags (with values > 255)
#define		MAESTRO_QUIT				256
#define		SCHEDULE_DATABASE			257
#define		SCHEDULE_BLOOM				258
#define		SCHEDULE_DOWNLOAD			259
#define		PROCESSOR_NAME				260

// When using the SRA "prefetch" command, we will log the error codes
// to track the failure modes
#define		PREFETCH_SUCCESS	0
#define		PREFETCH_NO_SRA		1000	// This is not a prefetch code, but is used to indicate no SRA file

struct search_by_accession
{
	inline bool operator()(const std::pair<SraAccession, size_t /*location*/> &m_lhs, const SraAccession &m_rhs) const
	{
		return m_lhs.first < m_rhs;
	};
};

struct BloomProgress
{
	#define BLOOM_PROGRESS_MEMBERS \
		VARIABLE(size_t, num_primary_align) \
		VARIABLE(size_t, curr_primary_align) \
		VARIABLE(size_t, num_unaligned_read) \
		VARIABLE(size_t, curr_unaligned_read) \
		VARIABLE(size_t, num_read) \
		VARIABLE(size_t, curr_read) \
		VARIABLE(size_t, curr_fragment) \
		VARIABLE(size_t, num_kmer) \
		VARIABLE(size_t, num_bp) \
		VARIABLE(size_t, log_2_counting_filter_len) \
		VARIABLE(std::string, error) \
		VARIABLE(bool, valid_read_collection) \

	#define VARIABLE(A, B) A B;
		BLOOM_PROGRESS_MEMBERS
	#undef VARIABLE

	BloomProgress()
	{
		// Since the error string cannot be assigned to zero, we need to manually
		// initialize all of the member variables here.
		num_primary_align = 0;
		curr_primary_align = 0;
		num_unaligned_read = 0;
		curr_unaligned_read = 0;
		num_read = 0;
		curr_read = 0;
		curr_fragment = 0;
		num_kmer = 0;
		num_bp = 0;
		log_2_counting_filter_len = 0;
		valid_read_collection = false;
	};

	template<class T> friend size_t mpi_size(const T &m_obj);
	template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, 
		T &m_obj);
	template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr,
		const T &m_obj);
};

template<> size_t mpi_size(const BloomProgress &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const BloomProgress &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, BloomProgress &m_obj);

// In maestro.cpp
void maestro_main(MaestroOptions &m_opt);
void worker_main(MaestroOptions &m_opt);

// In maestro_main.cpp
size_t estimate_num_bases(const std::string &m_path);

// In build_db.cpp
bool build_db(const std::string &m_filename, const BloomParam &m_param, const std::deque<std::string> &m_bloom_files);

// In make_bloom.cpp
// Return STATUS_BLOOM_SUCCESS --> For successful execution
// Return STATUS_BLOOM_FAIL --> For failure/error
// Return STATUS_BLOOM_INVALID --> Unable satisfy Bloom parameters
unsigned char make_bloom_filter(const SraAccession &m_acc, const FilterInfo &m_info, BloomParam &m_param,
	BloomProgress &m_progress, const std::string &m_bloom_dir, const MaestroOptions &m_opt,
	bool m_force_unaligned = false);

// In worder_main.cpp
bool remove_sra_files(const std::string &m_sra_dir, const std::string m_accession);

// In file_io.cpp
size_t get_accession_index(const SraAccession &m_accession, 
	const std::vector< std::pair<SraAccession, size_t /*location*/> > &m_accession_loc);
bool parse_accession_loc(std::vector< std::pair<SraAccession, size_t /*location*/> > &m_accession_loc, 
	const std::string &m_filename, bool m_verbose);
bool read_sra_repository(std::string &m_path);
bool write_status(const std::string &m_filename, unsigned char* m_status, const size_t &m_num_sra,
	size_t &m_database_index);
bool restore_status(const std::string &m_filename, unsigned char* m_status, const size_t &m_num_sra, 
	size_t &m_database_index, bool m_create_missing);

// In sra_meta.cpp
uint64_t number_of_bases(const std::string &m_accession);

#endif // __MAESTO