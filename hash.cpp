#include "hash.h"
#include <string.h>

using namespace std;

long murmur_hash(const string &m_seq, const unsigned int &m_seed);
long murmur_hash(const Word &m_seq, const unsigned int &m_kmer_len, 
	const unsigned int &m_seed);
	
inline unsigned int rotl32 ( unsigned int x, char r )
{
	return (x << r) | (x >> (32 - r));
}

inline unsigned int fmix ( unsigned int h )
{
	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;

	return h;
}

#define ROTL32(x,y) rotl32(x,y)

size_t bigsi_hash(const string &m_seq, const unsigned int &m_seed,
	const size_t &m_hash_range, const HashFunction &m_func)
{
	if(m_func == MURMUR_HASH){
		
		const long h = murmur_hash(m_seq, m_seed);
	
		//////////////////////////////////////////////////////////////////////////////
		// Note that Python handles the mod of negative numbers differently than C++!
		//////////////////////////////////////////////////////////////////////////////
		return ( (h < 0) ? ( (m_hash_range - 1) - (-h)%m_hash_range) : h%m_hash_range );
	}
	
	throw __FILE__ ":bigsi_hash: Unknown hash function";
	return 0;
}

size_t bigsi_hash(const Word &m_seq, const unsigned int &m_kmer_len,
	const unsigned int &m_seed, const size_t &m_hash_range, const HashFunction &m_func)
{
	if(m_func == MURMUR_HASH){
	
		const long h = murmur_hash(m_seq, m_kmer_len, m_seed);

		//////////////////////////////////////////////////////////////////////////////
		// Note that Python handles the mod of negative numbers differently than C++!
		//////////////////////////////////////////////////////////////////////////////
		return ( (h < 0) ? ( (m_hash_range - 1 ) - (-h)%m_hash_range) : h%m_hash_range );
	}
	
	throw __FILE__ ":bigsi_hash: Unknown hash function";
	return 0;
}

// From https://github.com/hajimes/mmh3/blob/master/MurmurHash3.cpp
long murmur_hash(const string &m_seq, const unsigned int &m_seed)
{
	const bool is_signed = true;
	
	const unsigned int block_size = 4;
	const unsigned int num_base = m_seq.size();
	const unsigned int nblocks = num_base/block_size;
		
	unsigned int h1 = m_seed;
	const unsigned int c1 = 0xcc9e2d51;
	const unsigned int c2 = 0x1b873593;
	
	// Body
	unsigned int offset = 0;
	
	for(unsigned int i = 0;i < nblocks;++i, offset += block_size){
		
		unsigned int k1 = (m_seq[offset + 3] << 24) |
				  (m_seq[offset + 2] << 16) |
				  (m_seq[offset + 1] << 8) |
				  (m_seq[offset + 0] << 0);
				  
		k1 *= c1;
		k1 = ROTL32(k1,15);
		k1 *= c2;

		h1 ^= k1;
		h1 = ROTL32(h1,13); 
		h1 = h1*5+0xe6546b64;
	}
	
	// tail
	unsigned int k1 = 0;
	
	switch(num_base & 3)
	{
		case 3: k1 ^= m_seq[offset + 2] << 16;
		case 2: k1 ^= m_seq[offset + 1] << 8;
		case 1: k1 ^= m_seq[offset];
		k1 *= c1; k1 = ROTL32(k1,15); k1 *= c2; h1 ^= k1;
	};

	//----------
	// finalization
	
	h1 ^= num_base;
	
	h1 = fmix(h1);
	
	int ret = 0;
	
	memcpy( &ret, &h1, sizeof(int) );
	
	return (is_signed ? 
		ret & 0xffffffffffffffff : 
		ret & 0x0ffffffff);
}

// From https://github.com/hajimes/mmh3/blob/master/MurmurHash3.cpp
long murmur_hash(const Word &m_seq, const unsigned int &m_kmer_len, 
	const unsigned int &m_seed)
{
	const bool is_signed = true;
	
	const unsigned int block_size = 4;
	const unsigned int nblocks = m_kmer_len/block_size;
		
	unsigned int h1 = m_seed;
	const unsigned int c1 = 0xcc9e2d51;
	const unsigned int c2 = 0x1b873593;
	
	// Body
	unsigned int offset = 0;
	
	#define BITS_TO_BASE(W,LEN,INDEX) \
		bits_to_base( ( W >> BITS_PER_BASE*( (LEN) - 1 - (INDEX) ) ) & 3 )
		
	for(unsigned int i = 0;i < nblocks;++i, offset += block_size){
		
		unsigned int k1 = (BITS_TO_BASE(m_seq, m_kmer_len, offset + 3) << 24) |
				  (BITS_TO_BASE(m_seq, m_kmer_len, offset + 2) << 16) |
				  (BITS_TO_BASE(m_seq, m_kmer_len, offset + 1) << 8) |
				  (BITS_TO_BASE(m_seq, m_kmer_len, offset + 0) << 0);
				  
		k1 *= c1;
		k1 = ROTL32(k1,15);
		k1 *= c2;

		h1 ^= k1;
		h1 = ROTL32(h1,13); 
		h1 = h1*5+0xe6546b64;
	}
	
	// tail
	unsigned int k1 = 0;
	
	switch(m_kmer_len & 3)
	{
		case 3: k1 ^= BITS_TO_BASE(m_seq, m_kmer_len, offset + 2) << 16;
		case 2: k1 ^= BITS_TO_BASE(m_seq, m_kmer_len, offset + 1) << 8;
		case 1: k1 ^= BITS_TO_BASE(m_seq, m_kmer_len, offset);
		k1 *= c1; k1 = ROTL32(k1,15); k1 *= c2; h1 ^= k1;
	};

	//----------
	// finalization
	
	h1 ^= m_kmer_len;
	
	h1 = fmix(h1);
	
	int ret = 0;
	
	memcpy( &ret, &h1, sizeof(int) );
	
	return (is_signed ? 
		ret & 0xffffffffffffffff : 
		ret & 0x0ffffffff);
}

string hash_name(const HashFunction &m_func)
{
	switch(m_func){
		case MURMUR_HASH:
			return std::string("murmur");
		case UNKNOWN_HASH:
			return std::string("unknown");
	};
	
	throw __FILE__ ":hash_name: Unknown hash function";
};

HashFunction parse_hash_function_name(const string &m_str)
{
	string name(m_str);
	
	for(string::iterator i = name.begin();i != name.end();++i){
		*i = tolower(*i);
	}
	
	for(int f = 0;f < UNKNOWN_HASH;++f){
	
		if( name == hash_name( HashFunction(f) ) ){
			return HashFunction(f);
		}
	}
	
	return UNKNOWN_HASH;
}
