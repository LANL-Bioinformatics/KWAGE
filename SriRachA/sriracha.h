#ifndef __SRIRACHA
#define __SRIRACHA

#include "mpi_util.h"
#include <unordered_map>

#define		MAP		unordered_map

#define		SRIRACHA_VERSION				"0.43"

#define		DEFAULT_MIN_READ_LENGTH			0
#define		MIN_KMER_LEN					3
#define		MAX_KMER_LEN					32

// The maximum number of matches to report for each SRA accession and query sequence
#define		DEFAULT_MAX_MATCH				100

// For the search-by-kmer strategy
#define		DEFAULT_KMER_LENGTH				11
#define		DEFAULT_KMER_MATCH_THRESHOLD	0.8
#define		DEFAULT_MIN_READ_COMPLEXITY		0.75
#define		DEFAULT_MIN_VALID_KMER			1

// The levels of verbosity
enum{
	SILENT,
	TACITERN,
	NORMAL,
	CHATTY
};

struct SearchMatch
{
	// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
	// ensure that structure variable are correctly serialized
	#define SEARCH_MATCH_MEMBERS \
		VARIABLE(unsigned int, read_index) \
		VARIABLE(unsigned int, read_subindex) \
		VARIABLE(float, score) \
		VARIABLE(std::string, read_seq)

	#define VARIABLE(A, B) A B;
		SEARCH_MATCH_MEMBERS
	#undef VARIABLE

	SearchMatch()
	{
		// Do nothing!
	};

	SearchMatch(const unsigned int &m_read_index, const unsigned int &m_read_subindex, const float &m_score,
		const std::string &m_seq):
		read_index(m_read_index), read_subindex(m_read_subindex), score(m_score), read_seq(m_seq)
	{

	};

	inline bool operator<(const SearchMatch &m_rhs) const
	{
		// If the scores are equal, sort matches by read index in
		// ascending order
		if(score == m_rhs.score){

			if(read_index == m_rhs.read_index){
				return (read_subindex < m_rhs.read_subindex);
			}

			return (read_index < m_rhs.read_index);
		}

		// Sort the matches by score in *descending* order
		return (score > m_rhs.score);
	};
};

template<> size_t mpi_size(const SearchMatch &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const SearchMatch &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, SearchMatch &m_obj);

// In search_by_kmer.cpp
void search_by_kmer(const std::string &m_seq, 
	const unsigned int &m_read_index, const unsigned int &m_read_subindex, void* m_param[]);

#endif // __SRIRACHA