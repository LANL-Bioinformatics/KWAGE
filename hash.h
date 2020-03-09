#ifndef __HASH_FUNCTION
#define __HASH_FUNCTION

#include <string>
#include "word.h"

enum {
	MURMUR_HASH,
	UNKNOWN_HASH // Must always be the last entry
};

typedef int HashFunction;

std::string hash_name(const HashFunction &m_func);
HashFunction parse_hash_function_name(const std::string &m_str);

size_t bigsi_hash(const std::string &m_seq, const unsigned int &m_seed,
	const size_t &m_hash_range, const HashFunction &m_func);

size_t bigsi_hash(const Word &m_seq, const unsigned int &m_kmer_len,
	const unsigned int &m_seed, const size_t &m_hash_range, 
	const HashFunction &m_func);
	
#endif // __HASH_FUNCTION
