#ifndef SRIRACHA_OPTIONS
#define SRIRACHA_OPTIONS

#include "mpi_util.h"

// We need to protect any commas that appear in template variables
// if we are going to combine them with X Macros
#define SINGLE_ARG(...) __VA_ARGS__

// The allowed search strategies
enum{
	SEARCH_BY_ALIGN,	// Seed and align
	SEARCH_BY_KMER,		// Kmer overlap
	SEARCH_BY_BLOOM,	// Bloom filter
	UNDEFINED_SEARCH
};

struct SrirachaOptions
{
	// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
	// ensure that structure variable are correctly serialized
	#define SRIRACHA_OPTION_MEMBERS \
		VARIABLE(std::deque<std::string>, sra_accession) \
		VARIABLE(std::deque<std::string>, input_sequence_files) \
		VARIABLE(std::string, output_filename) \
		VARIABLE(std::string, sra_accession_filename) \
		VARIABLE(unsigned int, kmer_len) \
		VARIABLE(float, kmer_match_threshold) \
		VARIABLE(float, min_read_complexity) \
		VARIABLE(unsigned int, min_read_length) \
		VARIABLE(unsigned int, min_valid_kmer) \
		VARIABLE(unsigned int, max_num_match) \
		VARIABLE(int, search_strategy) \
		VARIABLE(unsigned int, verbose) \
		VARIABLE(unsigned int, max_retry) \
		VARIABLE(bool, quit)

	#define VARIABLE(A, B) A B;
		SRIRACHA_OPTION_MEMBERS
	#undef VARIABLE

	SrirachaOptions()
	{
	};
	
	SrirachaOptions(int argc, char* argv[])
	{
		load(argc, argv);
	};
	
	void load(int argc, char* argv[]);
};

template<> size_t mpi_size(const SrirachaOptions &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const SrirachaOptions &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, SrirachaOptions &m_obj);

#endif // SRIRACHA_OPTIONS