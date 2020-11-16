#include "hash.h"
#include <string.h>
#include <immintrin.h>

using namespace std;

#define		NUM_SIMD_HASH	8

union SimdHash256 {
	
	SimdHash256()
	{
		// Do nothing
	};
	
	SimdHash256(const unsigned int &m_value)
	{
		simd = _mm256_set1_epi32(m_value);
	};

	~SimdHash256()
	{
		// Do nothing
	};
	
	__m256i simd;
	unsigned int v[NUM_SIMD_HASH];
};

const SimdHash256 five(5);

const unsigned int c1 = 0xcc9e2d51;
const unsigned int c2 = 0x1b873593;
const SimdHash256 c3(0xe6546b64);
const SimdHash256 c4(0x85ebca6b);
const SimdHash256 c5(0xc2b2ae35);

unsigned int murmur_hash32(const string &m_seq, const unsigned int &m_seed);
unsigned int murmur_hash32(const Word &m_seq, const unsigned int &m_kmer_len, 
	const unsigned int &m_seed);
void murmur_hash32(vector<size_t> &m_hash_values, const Word &m_seq, const unsigned int &m_kmer_len);

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

// The calling function must map the returned hash value to the required range
size_t bigsi_hash(const string &m_seq, const unsigned int &m_seed, const HashFunction &m_func)
{
	if(m_func == MURMUR_HASH_32){
				
		return size_t( murmur_hash32(m_seq, m_seed) );

		//////////////////////////////////////////////////////////////////////////////
		// Note that Python handles the mod of negative numbers differently than C++!
		//////////////////////////////////////////////////////////////////////////////
		//return ( (h < 0) ? ( (m_hash_range - 1) - (-h)%m_hash_range) : h%m_hash_range );
	}
	
	throw __FILE__ ":bigsi_hash: Unknown hash function";
	return 0;
}

// The calling function must map the returned hash value to the required range
size_t bigsi_hash(const Word &m_seq, const unsigned int &m_kmer_len,
	const unsigned int &m_seed, const HashFunction &m_func)
{
	if(m_func == MURMUR_HASH_32){
	
		return size_t( murmur_hash32(m_seq, m_kmer_len, m_seed) );

		//////////////////////////////////////////////////////////////////////////////
		// Note that Python handles the mod of negative numbers differently than C++!
		//////////////////////////////////////////////////////////////////////////////
		//return ( (h < 0) ? ( (m_hash_range - 1 ) - (-h)%m_hash_range) : h%m_hash_range );
	}
	
	throw __FILE__ ":bigsi_hash: Unknown hash function";
	return 0;
}

// The calling function must map the returned hash value to the required range
void bigsi_hash(vector<size_t> &m_hash_values, const Word &m_seq, const unsigned int &m_kmer_len,
	const HashFunction &m_func)
{
	if(m_func == MURMUR_HASH_32){
	
		murmur_hash32(m_hash_values, m_seq, m_kmer_len);

		return;
	}
	
	throw __FILE__ ":bigsi_hash: Unknown hash function";
}

// This is the 32-bit Murmur3 hash function from:
// https://github.com/hajimes/mmh3/blob/master/MurmurHash3.cpp
// This function is included as a reference, but is not actually called. Instead,
// we used a SIMD-parallel version of the MurmurHash3 function (see below).
unsigned int murmur_hash32(const string &m_seq, const unsigned int &m_seed)
{	
	const unsigned int block_size = 4;
	const unsigned int num_base = m_seq.size();
	const unsigned int nblocks = num_base/block_size;
		
	unsigned int h1 = m_seed;
	//const unsigned int c1 = 0xcc9e2d51;
	//const unsigned int c2 = 0x1b873593;
	
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
	
	//int ret = 0;
	
	//memcpy( &ret, &h1, sizeof(int) );
	
	//return (IS_SIGNED ? 
	//	ret & 0xffffffffffffffff : 
	//	ret & 0x0ffffffff);

	return h1;
}

// This is the 32-bit Murmur3 hash function from:
// https://github.com/hajimes/mmh3/blob/master/MurmurHash3.cpp
// This function is included as a reference, but is not actually called. Instead,
// we used a SIMD-parallel version of the MurmurHash3 function (see below).
unsigned int murmur_hash32(const Word &m_seq, const unsigned int &m_kmer_len, 
	const unsigned int &m_seed)
{	
	const unsigned int block_size = 4;
	const unsigned int nblocks = m_kmer_len/block_size;
		
	unsigned int h1 = m_seed;
	//const unsigned int c1 = 0xcc9e2d51;
	//const unsigned int c2 = 0x1b873593;
	
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
	
	//int ret = 0;
	
	//memcpy( &ret, &h1, sizeof(int) );
	
	//return (IS_SIGNED ? 
	//	ret & 0xffffffffffffffff : 
	//	ret & 0x0ffffffff);
	return h1;
}

// This is based on the 32-bit Murmur3 hash function from:
// https://github.com/hajimes/mmh3/blob/master/MurmurHash3.cpp
// SIMD instructions are used to compute multiple hash values in parallel.
void murmur_hash32(vector<size_t> &m_hash_values, const Word &m_seq, const unsigned int &m_kmer_len)
{
	const unsigned int num_seed = m_hash_values.size();

	if(num_seed > NUM_SIMD_HASH){
		throw __FILE__ ":murmur_hash32: Exceeded maximum number of hash values";
	}

	const unsigned int block_size = 4;
	const unsigned int nblocks = m_kmer_len/block_size;
	
	SimdHash256 h1;

	for(unsigned int seed = 0;seed < num_seed;++seed){
		h1.v[seed] = seed;
	}

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

		// h1 ^= k1
		h1.simd = _mm256_xor_si256( h1.simd, _mm256_set1_epi32(k1) );

		// h1 = ROTL32(h1,13); --> return (h1 << r) | (h1 >> (32 - r));
		//h1.simd = _mm256_rol_epi32(h1.simd, 13); // <-- requires avx512vl (not sufficiently portable yet)
		h1.simd = _mm256_or_si256( _mm256_slli_epi32(h1.simd, 13), _mm256_srli_epi32(h1.simd, 32 - 13) );

		// h1 = h1*5 + 0xe6546b64;
		h1.simd = _mm256_mullo_epi32(h1.simd, five.simd); // h1 += h1*5
		h1.simd = _mm256_add_epi32(h1.simd, c3.simd); // h1 += 0xe6546b64	
	}

	// tail
	unsigned int k1 = 0;
	
	switch(m_kmer_len & 3)
	{
		case 3: k1 ^= BITS_TO_BASE(m_seq, m_kmer_len, offset + 2) << 16;
		case 2: k1 ^= BITS_TO_BASE(m_seq, m_kmer_len, offset + 1) << 8;
		case 1: k1 ^= BITS_TO_BASE(m_seq, m_kmer_len, offset);
		k1 *= c1; k1 = ROTL32(k1,15); k1 *= c2;

		// h1 ^= k1
		h1.simd = _mm256_xor_si256( h1.simd, _mm256_set1_epi32(k1) );
	};

	//----------
	// finalization
	
	// h1 ^= m_kmer_len;
	h1.simd = _mm256_xor_si256( h1.simd, _mm256_set1_epi32(m_kmer_len) );

	// h1 = fmix(h1), where fmix is:
	//h ^= h >> 16;
	h1.simd = _mm256_xor_si256( h1.simd, _mm256_srli_epi32(h1.simd, 16) );

	//h *= 0x85ebca6b;
	h1.simd = _mm256_mullo_epi32(h1.simd, c4.simd);

	//h ^= h >> 13;
	h1.simd = _mm256_xor_si256( h1.simd, _mm256_srli_epi32(h1.simd, 13) );

	//h *= 0xc2b2ae35;
	h1.simd = _mm256_mullo_epi32(h1.simd, c5.simd);

	//h ^= h >> 16;
	h1.simd = _mm256_xor_si256( h1.simd, _mm256_srli_epi32(h1.simd, 16) );

	// While there is no AVX SIMD modulo division operator, can we build one
	// via: X%Y = X - (X/Y)*Y ??
	// Unfortunately, there does not appear to be an intrinsic for
	// integer division. While _mm256_div_epi32() is mentioned on the
	// Intel website:https://software.intel.com/sites/landingpage/IntrinsicsGuide
	// the clang compiler doesn't support it!

	for(unsigned int seed = 0;seed < num_seed;++seed){
		m_hash_values[seed] = h1.v[seed];
	}
}

string hash_name(const HashFunction &m_func)
{
	switch(m_func){
		case MURMUR_HASH_32:
			return std::string("murmur32");
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
