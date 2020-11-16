#ifndef __HASH_FUNCTION
#define __HASH_FUNCTION

#include <string>
#include <vector>
#include "word.h"

enum {
	MURMUR_HASH_32, // <-- Make it clear that this is a 32bit hash function!
	UNKNOWN_HASH // Must always be the last entry
};

typedef int HashFunction;

std::string hash_name(const HashFunction &m_func);
HashFunction parse_hash_function_name(const std::string &m_str);

// Currently, we only support 32 bit hash functions, which limits the 
// Bloom filter size to a length of 2^32 bits
size_t bigsi_hash(const std::string &m_seq, const unsigned int &m_seed,
	const HashFunction &m_func);

size_t bigsi_hash(const Word &m_seq, const unsigned int &m_kmer_len,
	const unsigned int &m_seed, 
	const HashFunction &m_func);

// Compute the vector of hash values, where the same value (i.e. m_seq) is hashed
// will a different seed that is equal to the index of the m_hash_values vector
void bigsi_hash(std::vector<size_t> &m_hash_values, const Word &m_seq, const unsigned int &m_kmer_len,
	const HashFunction &m_func);

#endif // __HASH_FUNCTION
