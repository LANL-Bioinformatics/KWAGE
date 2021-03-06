#include <iostream>
#include "word.h"
#include "hash.h"
#include "bloom.h"
#include <ncbi-vdb/NGS.hpp> // For openReadCollection

// For non-temporal memory access that bypasses the cache
#include <smmintrin.h>
#include <emmintrin.h>

using namespace std;

#define		MAX_NUM_HASH	5

#define		LOG_COUNT_FILTER_LEN	32

#define		NUM_COUNT_HASH			2
#define		MAX_COUNT				15

#define COUNT_MAP	std::unordered_map

// Pack two counting Bloom filters into 8 bits.
// Each has a max count of 2^4 - 1 = 15
struct CountingBloom
{
	typedef unsigned char CountingBloomBaseType;

	CountingBloomBaseType first : 4;
	CountingBloomBaseType second : 4;
};

void count_words(const COUNT_MAP<Word, unsigned int> &m_buffer, const unsigned int &m_min_kmer_count,
	const unsigned int &m_kmer_len, const HashFunction &m_hash_func, CountingBloom *m_count_ptr,
	size_t &m_num_valid_kmer, vector<BitVector> &m_unique_bits, const size_t &m_hash_count_mask, 
	const size_t &m_hash_seq_mask);

int main(int argc, char *argv[])
{
	try{
		if(argc != 2){

			cerr << "Usage: " << argv[0] << " <SRA accession>" << endl;
			return EXIT_SUCCESS;
		}

		time_t profile = time(NULL);
		const time_t start_profile = profile;

		const string acc = argv[1];

		const unsigned int min_log_2_filter_len = 18;
		const unsigned int max_log_2_filter_len = 32;
		const unsigned int kmer_len = 31;
		const HashFunction hash_func = MURMUR_HASH_32;
		const float false_positive_probability = 0.25;

		const unsigned int min_kmer_count = 5;

		//const size_t max_num_read = 1000000;
		const size_t max_num_read = 2000000;

		CountingBloom *bcount = NULL;

		// DEBUG
		cerr << "sizeof(CountingBloom) = " << sizeof(CountingBloom) << endl;
		cerr << "NUM_COUNT_HASH = " << NUM_COUNT_HASH << endl;
		cerr << "MAX_COUNT = " << MAX_COUNT << endl;

		cerr << "approximate_max_kmers = " << 
			approximate_max_kmers(false_positive_probability,
				hash_func, min_log_2_filter_len, max_log_2_filter_len) << endl;

		// The size of the counting bloom filter is the size of the maximum allowed Bloom filter
		const size_t num_seq_bloom = 1ULL << max_log_2_filter_len;
		const size_t num_count_bloom = 1ULL << LOG_COUNT_FILTER_LEN;

		// Since we are restricting the Bloom filter lengths to be a power of
		// two, we can use the following replacement for modulo division:
		// X % (2^n) = X & (2^n - 1)

		size_t hash_seq_mask = 0;
		size_t hash_count_mask = 0;

		for(size_t i = 0;i < max_log_2_filter_len;++i){
			hash_seq_mask |= (1ULL << i);
		}

		for(size_t i = 0;i < LOG_COUNT_FILTER_LEN;++i){
			hash_count_mask |= (1ULL << i);
		}

		size_t num_valid_kmer = 0;

		//bcount = new CountingBloom[num_count_bloom];
		posix_memalign( (void**)&bcount, 16, num_count_bloom*sizeof(CountingBloom) );

		if(bcount == NULL){
			throw __FILE__ ":make_bloom_filter: Unable to allocate the counting bloom filter";
		}

		cerr << "Allocated (" << std::hex << size_t(bcount) << std::dec <<  ")" << endl;

		// The counting Bloom filter must be initialized to zero before we add any kmers
		memset(bcount, 0, num_count_bloom);

		vector<BitVector> unique_bits( MAX_NUM_HASH, BitVector(num_seq_bloom) );

		for(size_t h = 0;h < MAX_NUM_HASH;++h){
			unique_bits[h].unset_all_bits();
		}

		#define	GROUND_TRUTH

		#ifdef GROUND_TRUTH
		vector< deque<Word> > curr_kmers(1);
		#endif // GROUND_TRUTH

		// Digest the input sequence into kmers and insert each kmer into the counting 
		// Bloom filter
		ngs::ReadCollection run(  ncbi::NGS::openReadCollection(acc) );
			
		// Note that num_read is the number of either paired or
		// unpaired reads. For paired reads, this is half the
		// the number of sequences!
		const size_t num_read = run.getReadCount(ngs::Read::all);

		cerr << "Found " << num_read << " reads" << endl;

		// The ReadRange is 1's based
		ngs::ReadIterator run_iter = 
			ngs::ReadIterator( run.getReadRange ( 1, min(max_num_read, num_read), ngs::Read::all ) );
			
		size_t seq_count = 0;
		size_t read_count = 0;

		cerr << "Reading ..." << endl;

		double begin_read_time = MPI_Wtime();

		// Count words in batches to reduce the number of times we hash and count the *same*
		// kmer multiple times.
		// - SRA records are sorted for good compression, which increases the chances that
		//	 adjacent reads will have similar kmer compositions.
		// - Hashing and Bloom filter memory access is expensive, so count kmers in batches
		//   that are small enough to store in a Map.
		// - Kmers that appear multiple times will only be hashed and counted once.
		//const size_t max_word_buffer = 10000; // 136 sec
		//const size_t max_word_buffer = 5000; // 130 sec
		const size_t max_word_buffer = 2000; // 129 sec
		//const size_t max_word_buffer = 1000; //  132 sec
	
		COUNT_MAP<Word, unsigned int> word_buffer;

		while( run_iter.nextRead() ){

			++read_count;

			if(word_buffer.size() >= max_word_buffer){

				count_words(word_buffer, min_kmer_count, kmer_len, hash_func, bcount, num_valid_kmer, 
					unique_bits, hash_count_mask, hash_seq_mask);

				word_buffer.clear();
			}
			
			while( run_iter.nextFragment() ){
				
				++seq_count;
										
				const ngs::StringRef &seq = run_iter.getFragmentBases();

				ForEachDuplexWord(seq.data(), seq.data() + seq.size(), kmer_len)

					if(ValidWord){

						const Word w = CanonicalWord;

						#ifdef GROUND_TRUTH
						curr_kmers[0].push_back(w);
						#endif // GROUND_TRUTH

						++word_buffer[w];
						//debug.insert(w);
						//++debug_num_word;
					}
				EndWord
			}

			// DEBUG
			const size_t update_every = 100000;

			if(read_count%update_every == 0){

				const double end_read_time = MPI_Wtime();

				//cerr << end_read_time - begin_read_time << " sec; " << (100.0*debug.size())/debug_num_word << '%' << endl;
				cerr << end_read_time - begin_read_time << " sec" << endl;

				begin_read_time = end_read_time;
			}
		}

		count_words(word_buffer, min_kmer_count, kmer_len, hash_func, bcount, num_valid_kmer, 
			unique_bits, hash_count_mask, hash_seq_mask);

		word_buffer.clear();

		deque<Word> valid_kmers; // kmers that have passed the minimum count threshold
		
		#ifdef GROUND_TRUTH
		// Count the occurance of each k-mer for frequency-based
		// k-mer filtering (to remove sequencing errors)
		find_abundant_kmers(valid_kmers, curr_kmers, min_kmer_count);

		// Discard any kmers remaining in curr_kmers, they did not appear
		// at least opt.min_kmer_count times
		curr_kmers.clear(); // Free memory
		#endif // GROUND_TRUTH

		// DEBUG
		cerr << "Counted kmers in " << time(NULL) - profile << " sec" << endl;
		cerr << "brute force num unique kmer = " << valid_kmers.size() << endl;
		cerr << "num_valid_kmer = " << num_valid_kmer << endl;

		const size_t max_unique_kmer = max(valid_kmers.size(), num_valid_kmer);

		BloomParam param;

		try{

			param = optimal_bloom_param(kmer_len,
				max_unique_kmer, 
				false_positive_probability,
				hash_func, 
				min_log_2_filter_len, 
				max_log_2_filter_len);

			cerr << "param.log_2_filter_len = " << param.log_2_filter_len << endl;
			cerr << "param.num_hash = " << param.num_hash << endl;
		}
		catch(...){

			// We were unable to find Bloom filter parameters that satisfied the 
			//requested false_positive_probability.

			if(bcount != NULL){

				//delete [] bcount;
				free(bcount);
				bcount = NULL;
			}

			throw __FILE__ ": Unable to find valid Bloom filter parameters";
		}

		// The filter_len (i.e. number of bits in a Bloom filter) can be very large.
		// Use a 64-bit unsigned integer for now ...
		const size_t filter_len = param.filter_len();
		
		cerr << "filter len = " << filter_len << endl;

		BloomFilter filter(param);
		BloomFilter raw_filter(param);

		filter.unset_all_bits();
		raw_filter.unset_all_bits();

		for(deque<Word>::const_iterator i = valid_kmers.begin();i != valid_kmers.end();++i){

			// The BIGSI python implementation starts from 0 when seeding the
			// hash function
			for(size_t h = 0;h < param.num_hash;++h){
				raw_filter.set_bit( bigsi_hash(*i, kmer_len, h, hash_func)%filter_len );
			}
		}

		if(bcount != NULL){

			profile = time(NULL);

			BitVector::BLOCK *dst_ptr = filter.ptr();
			const uint64_t num_dst_block = filter.num_block();

			for(size_t h = 0;h < param.num_hash;++h){

				BitVector::BLOCK *src_ptr = unique_bits[h].ptr();
				const uint64_t num_src_block = unique_bits[h].num_block();

				for(uint64_t i = 0;i < num_src_block;i += num_dst_block){

					for(uint64_t j = 0;j < num_dst_block;++j){
						dst_ptr[j] |= src_ptr[i + j];
					}
				}
			}

			profile = time(NULL) - profile;

			cerr << "Set bits in the final filter in " << profile << " sec" << endl;
		}

		// Clean up the counting filter
		if(bcount != NULL){
			
			delete [] bcount;
			bcount = NULL;
		}	

		size_t num_discord = 0;

		for(size_t i = 0;i < filter_len;++i){

			if( filter.get_bit(i) != raw_filter.get_bit(i) ){

				++num_discord;

				// Tests show that collapsing the CountingBloom filter always sets *extra* bits
				// (as expected).
				//cout << i << '\t' << int( filter.get_bit(i) ) << '\t' <<raw_filter.get_bit(i) << endl;
			}
		}

		cerr << num_discord << " of " << filter_len << " bits (" 
			<< (100.0*num_discord)/filter_len << "%) disagree" << endl;

		profile = time(NULL) - start_profile;
		
		cerr << "Tested counting Bloom filter in " << profile << " sec" << endl;

	}
	catch(const char* error){

		cerr << "Caught the error: " << error << endl;
		return EXIT_FAILURE;
	}
	catch(...){

		cerr << "Caught an unhandled error" << endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

#define		MEM_ALIGN		16ULL

union Buffer128
{
	__m128i simd;
	unsigned char byte[16];
} __attribute((aligned(MEM_ALIGN)));

void count_words(const COUNT_MAP<Word, unsigned int> &m_buffer, const unsigned int &m_min_kmer_count,
	const unsigned int &m_kmer_len, const HashFunction &m_hash_func, CountingBloom *m_count_ptr, 
	size_t &m_num_valid_kmer, vector<BitVector> &m_unique_bits, const size_t &m_hash_count_mask, 
	const size_t &m_hash_seq_mask)
{

	// The counting Bloom filter and the sequence Bloom filter can have different
	// ranges, so we will map the hash values to the required range as needed below.
	//size_t total = 0;

	vector<size_t> hash_values(MAX_NUM_HASH);

	// No speed up observed for the CACHE_BYPASS. Off by a few bits as well (not sure why!)
	//#define	CACHE_BYPASS

	#ifdef CACHE_BYPASS
	Buffer128 local;
	#endif // CACHE_BYPASS

	for(COUNT_MAP<Word, unsigned int>::const_iterator kmer_iter = m_buffer.begin();kmer_iter != m_buffer.end();++kmer_iter){

		//total += kmer_iter->second;

		bigsi_hash(hash_values, kmer_iter->first, m_kmer_len, m_hash_func);

		#ifdef CACHE_BYPASS
		unsigned char min_count = MAX_COUNT;

		for(int i = 0;i < NUM_COUNT_HASH;++i){

			size_t hash_index = hash_values[i] & m_hash_count_mask;
			size_t local_index = hash_index%MEM_ALIGN;

			__m128i* ptr = (__m128i*)(m_count_ptr + hash_index - hash_index%MEM_ALIGN);

			local.simd = _mm_stream_load_si128(ptr);

			min_count = min(min_count, local.byte[local_index]);

			if(local.byte[local_index] < m_min_kmer_count){

				local.byte[local_index] = min( MAX_COUNT, int(local.byte[local_index] + kmer_iter->second) );
				_mm_stream_si128(ptr, local.simd);
			}
		}
		#endif // CACHE_BYPASS

		#define ORIGINAL_VERSION
		#ifdef ORIGINAL_VERSION

		//if( (hash_values[0] & m_hash_count_mask) >= ( 1ULL << LOG_COUNT_FILTER_LEN) ){
		//	throw __FILE__ ":count_words: Overflow 0";
		//}

		//if( (hash_values[1] & m_hash_count_mask) >= ( 1ULL << LOG_COUNT_FILTER_LEN) ){
		//	throw __FILE__ ":count_words: Overflow 1";
		//}

		__builtin_prefetch(m_count_ptr + (hash_values[1] & m_hash_count_mask), 0, 0);

		CountingBloom *ptr = m_count_ptr + (hash_values[0] & m_hash_count_mask);
		
		unsigned char min_count = ptr->first;

		if(ptr->first < MAX_COUNT){
			ptr->first = min( MAX_COUNT, int(ptr->first + kmer_iter->second) );
		}

		ptr = m_count_ptr + (hash_values[1] & m_hash_count_mask);
		
		min_count = min(min_count, ptr->second);

		if(ptr->second < MAX_COUNT){
			ptr->second = min( MAX_COUNT, int(ptr->second + kmer_iter->second) );
		}

		#endif // ORIGINAL_VERSION

		// The count is clamped to be <= MAX_COUNT. Count a kmer as "valid" when it's minimum
		// count crosses the the m_min_kmer_count threshold
		if( (min_count < m_min_kmer_count) && ( (min_count + kmer_iter->second) >= m_min_kmer_count ) ){

			++m_num_valid_kmer;

			for(size_t h = 0;h < MAX_NUM_HASH;++h){
				m_unique_bits[h].set_bit(hash_values[h] & m_hash_seq_mask);
			}
		}
	}

	//cout << (100.0*m_buffer.size())/total << endl;
}