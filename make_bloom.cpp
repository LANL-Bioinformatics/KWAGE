// Convert SRA sequence files and metadata into Bloom filters
#include <iostream>

#include <math.h>

#include <ncbi-vdb/NGS.hpp> // For openReadCollection

#include "maestro.h"
#include "word.h"
#include "file_util.h"
#include "mem_usage.h"

using namespace std;

extern int mpi_rank; // for debugging

// The size of the counting Bloom filter, 2^LOG_COUNT_FILTER_LEN, is independent of the
// length of the Bloom filters for storing SRA data. The counting filter can consume
// a *lot* of RAM -- exactly 2^LOG_COUNT_FILTER_LEN * sizeof(CountingBloom) bytes.
//
#define		MAX_LOG_COUNT_FILTER_LEN	32ULL	// A value greater than 32 will require a 64-bit hash function
#define		MIN_LOG_COUNT_FILTER_LEN	18ULL	

// The maximum allowed false positive probability for the counting Bloom filters
#define		COUNT_FILTER_FP				1.0e-2

// Use a total of NUM_COUNT_HASH hash functions when constructing counting Bloom filters to count kmers
// and identify the "valid" kmers whose count is greater than, or equal to, the
// user-supplied threshold.
//
// Two counting Bloom filters are used. Each use two hash functions (our goal is to
// minimize the false positive rate as a function of the *total* number of hash function evaluations).
// Two and two is predicted to have slightly better performance than one and three for all but the largest
// filters
// Here are benchmarks for one and 'N' approach:
//		- The first is a single hash function filter (for the worst case scenario of 
//		  the number of kmers >= 2^32)
//		- The second is a 5 hash function filter (for the scenarios when the number of
//		   kmers <= 2^29)
//			- For DRR000347, Bloom log_2_filter_len = 30, num hash = 4
//			- Using 5 hash function for filter 2 --> 1481 sec, 0.111901% bit difference to 2xFive Hash filter
//			- Using 4 hash function for filter 2 --> 1417 sec, 0.0990211% bit difference to 2xFive Hash filter
//			- Using 3 hash function for filter 2 --> 1300 sec, 0.113448% bit difference to 2xFive Hash filter
//			- Using 2 hash function for filter 2 --> 1238 sec, 0.132031% bit difference to 2xFive Hash filter
//			- Using 1 hash function for filter 2 --> 1141 sec, 0.221165% bit difference to 2xFive Hash filter

// Single hash function counting filters are used because they are optimal in
// the *worst* case scenario of a very large number of kmers. They are suboptimal
// for smaller numbers of kmers.
#define		NUM_COUNT_HASH			4 // 2 for the first filter + 2 for the second filter

#if (MAX_NUM_HASH > 8)
#error The current bigsi++ hash function supports a maximum of 8 hash values/word
#endif

// If we need to count higher to determine valid kmers, we will need more bits! 
// We can get more bits for counting by:
//  a) Increase the number of available bits by changing CountingBloomBaseType to be
//	   unsigned short (or the type of you choice), or
//  b) Reducing the number of counting Bloom filters to one
#define		MAX_COUNT				15U

struct CountingBloom
{
	typedef unsigned char CountingBloomBaseType;

	CountingBloomBaseType first : 4;	// Store the count of the first, single hash filter
	CountingBloomBaseType second : 4;	// Store the count of the second, single hash filter
};

void count_words(CountingBloom *m_count_ptr, vector<BitVector> &m_valid_bits, 
	size_t &m_num_valid_kmer, const ngs::StringRef &m_seq, 
	const size_t &m_hash_seq_mask, const size_t &m_hash_count_mask, 
	const MaestroOptions &m_opt);

unsigned char make_bloom_filter(const SraAccession &m_acc, const FilterInfo &m_info, BloomParam &m_param,
	BloomProgress &m_progress, 
	const string &m_bloom_dir, const MaestroOptions &m_opt, bool m_force_unaligned /*= false*/)
{
	//#define DEBUG_BLOOM

	#ifdef DEBUG_BLOOM
	cerr << "[" << mpi_rank << "] in make_bloom_filter" << endl;
	#endif // DEBUG_BLOOM

	CountingBloom *bcount = NULL;

	try{

		if(m_opt.min_kmer_count > MAX_COUNT){
			throw __FILE__ ":make_bloom_filter: min_kmer_count is too large. See the comments in make_bloom.cpp for parameter settings.";
		}

		if( (MAX_LOG_COUNT_FILTER_LEN > 32) && (m_opt.hash_func == MURMUR_HASH_32) ){
			throw __FILE__ ":make_bloom_filter: The MAX_LOG_COUNT_FILTER_LEN is too large for a 32-bit hash function.";
		}

		assert(MAX_LOG_COUNT_FILTER_LEN > MIN_LOG_COUNT_FILTER_LEN);

		const string accession = accession_to_str(m_acc);

		// The sizes the the counting Bloom filter and the sequence Bloom filter are independent
		const size_t max_seq_bloom = 1ULL << m_opt.max_log_2_filter_len;

		// Assume the worse-case scenario and use the largest allowed counting Bloom filter(s)
		m_progress.log_2_counting_filter_len = MAX_LOG_COUNT_FILTER_LEN;

		// Try to get the number of base pairs from the SRA metadata
		const uint64_t num_bp = number_of_bases(accession);

		// Did we obtain a valid number of base pairs from the SRA metadata?
		if(num_bp > 0){

			// Required Counting Bloom filter length assuming two counting filters, each with
			// two hash functions (and the same number of bits).
			const double counting_length = 1.0/( 1.0 - pow( 1.0 - pow(COUNT_FILTER_FP, 1.0/4.0) , 1.0/(2*num_bp) ) );

			m_progress.log_2_counting_filter_len = ceil( log( counting_length)/log(2.0) );

			// Clamp the log_count_filter_len to be within the 
			// allowed range: [MIN_LOG_COUNT_FILTER_LEN, MAX_LOG_COUNT_FILTER_LEN]
			if(m_progress.log_2_counting_filter_len > MAX_LOG_COUNT_FILTER_LEN){
				m_progress.log_2_counting_filter_len = MAX_LOG_COUNT_FILTER_LEN;
			}

			if(m_progress.log_2_counting_filter_len < MIN_LOG_COUNT_FILTER_LEN){
				m_progress.log_2_counting_filter_len = MIN_LOG_COUNT_FILTER_LEN;
			}
		}

		const size_t num_count_bloom = 1ULL << m_progress.log_2_counting_filter_len;

		// Since we are restricting the Bloom filter lengths to be a power of
		// two, we can use the following replacement for modulo division:
		// X % (2^n) = X & (2^n - 1)

		size_t hash_seq_mask = 0;
		size_t hash_count_mask = 0;

		for(size_t i = 0;i < m_opt.max_log_2_filter_len;++i){
			hash_seq_mask |= (1ULL << i);
		}

		for(size_t i = 0;i < m_progress.log_2_counting_filter_len;++i){
			hash_count_mask |= (1ULL << i);
		}

		const size_t max_num_kmer = approximate_max_kmers(m_opt.false_positive_probability,
			m_opt.hash_func, m_opt.min_log_2_filter_len, m_opt.max_log_2_filter_len);

		bcount = new CountingBloom[num_count_bloom];

		if(bcount == NULL){
			throw __FILE__ ":make_bloom_filter: Unable to allocate the counting bloom filter";
		}

		// The counting Bloom filter must be initialized to zero before we count any kmers
		memset(bcount, 0, num_count_bloom);

		// Accumulate the per-hash bits for each kmer that appears at
		// least m_opt.min_kmer_count times.
		vector<BitVector> valid_bits( MAX_NUM_HASH, BitVector(max_seq_bloom) );

		for(size_t h = 0;h < MAX_NUM_HASH;++h){
			valid_bits[h].unset_all_bits();
		}

		// Digest the input sequence into kmers and insert each kmer into the counting 
		// Bloom filter
		ngs::ReadCollection run(  ncbi::NGS::openReadCollection(accession) );
		
		m_progress.valid_read_collection = true;

		// We are following the advice of Kurt Rodamer and Kenneth Durbrow @ NCBI
		// 		** This approach will *miss* a small number of reads that are only
		// 		** partially aligned (i.e. only one read of a pair aligned). If all reads
		// 		** in an SRA record are needed, then iterate through all reads using
		// 		** ngs::ReadIterator( run.getReadRange ( 1, num_read, ngs::Read::all ) );
		//
		// Step 1: Does the SRA record contain aligned reads?
		m_progress.num_primary_align = run.getAlignmentCount(ngs::Alignment::primaryAlignment);

		#ifdef DEBUG_BLOOM
		cerr << "[" << mpi_rank << "] found " << m_progress.num_primary_align << " primary alignments in " 
			<< accession_to_str(m_acc) << endl;
		#endif // DEBUG_BLOOM

		if( (m_progress.num_primary_align > 0) && !m_force_unaligned){

			// Step 2: Read the primaryAlignment sequences first, as this minimized the amount
			// of random I/O. If needed, we can also start from a particular alignment using:
			// getAlignmentRange(uint64_t first, uint64_t count, Alignment::AlignmentCategory categories)
			// See ngs/ngs-sdk/include/ngs/ReadCollection.hpp for more member functions.
			ngs::AlignmentIterator align_iter = run.getAlignments(ngs::Alignment::primaryAlignment);

			while(align_iter.nextAlignment() ){

				// Track the total number of bases read
				m_progress.num_bp += align_iter.getAlignedFragmentBases().size();

				count_words(bcount, valid_bits, m_progress.num_kmer, 
					align_iter.getAlignedFragmentBases(), 
					hash_seq_mask, hash_count_mask, m_opt);

				// Track the progress for restarting and/or error reporting
				++m_progress.curr_primary_align;

				if(max_num_kmer < m_progress.num_kmer){

					delete [] bcount;
					bcount = NULL;
					
					return STATUS_BLOOM_INVALID;
				}
			}

			// Step 3: Read the unaligned sequences
			m_progress.num_unaligned_read = run.getReadCount(ngs::Read::unaligned);

			#ifdef DEBUG_BLOOM
			cerr << "[" << mpi_rank << "] found " << m_progress.num_unaligned_read << " unaligned reads in " 
				<< accession_to_str(m_acc) << endl;
			#endif // DEBUG_BLOOM

			if(m_progress.num_unaligned_read > 0){

				// Need to use getReads() -- getReadRange() does not appear to work for unaligned reads
				ngs::ReadIterator read_iter = ngs::ReadIterator( run.getReads(ngs::Read::unaligned) );

				while( read_iter.nextRead() ){
					
					m_progress.curr_fragment = 0;

					while( read_iter.nextFragment() ){

						// Track the total number of bases read
						m_progress.num_bp += read_iter.getFragmentBases().size();

						count_words(bcount, valid_bits, m_progress.num_kmer, 
							read_iter.getFragmentBases(), 
							hash_seq_mask, hash_count_mask, m_opt);
						
						// Track the progress for restarting and/or error reporting
						++m_progress.curr_fragment;

						if(max_num_kmer < m_progress.num_kmer){

							delete [] bcount;
							bcount = NULL;
							
							return STATUS_BLOOM_INVALID;
						}
					}

					// Track the progress for restarting and/or error reporting
					++m_progress.curr_unaligned_read;
				}				
			}
		}
		else{ // We not *not* have aligned reads

			m_progress.num_read = run.getReadCount(ngs::Read::all);

			#ifdef DEBUG_BLOOM
			cerr << "[" << mpi_rank << "] found " << m_progress.num_read << " reads in " 
				<< accession_to_str(m_acc) << endl;
			#endif // DEBUG_BLOOM

			ngs::ReadIterator read_iter = 
				ngs::ReadIterator( run.getReadRange ( 1, m_progress.num_read, ngs::Read::all ) );

			while( read_iter.nextRead() ){
				
				m_progress.curr_fragment = 0;

				while( read_iter.nextFragment() ){
					
					// Track the total number of bases read
					m_progress.num_bp += read_iter.getFragmentBases().size();

					count_words(bcount, valid_bits, m_progress.num_kmer, 
						read_iter.getFragmentBases(), 
						hash_seq_mask, hash_count_mask, m_opt);
					
					// Track the progress for restarting and/or error reporting
					++m_progress.curr_fragment;

					if(max_num_kmer < m_progress.num_kmer){

						delete [] bcount;
						bcount = NULL;
						
						return STATUS_BLOOM_INVALID;
					}
				}

				// Track the progress for restarting and/or error reporting
				++m_progress.curr_read;
			}
		}

		#ifdef DEBUG_BLOOM
		cerr << "[" << mpi_rank << "] finished digesting reads" << endl;
		#endif // DEBUG_BLOOM
			
		try{

			m_param = optimal_bloom_param(m_opt.kmer_len,
				m_progress.num_kmer, 
				m_opt.false_positive_probability,
				m_opt.hash_func, 
				m_opt.min_log_2_filter_len, 
				m_opt.max_log_2_filter_len);
			
			#ifdef DEBUG_BLOOM
			cerr << "[" << mpi_rank << "] found valid Bloom parameters" << endl;
			#endif // DEBUG_BLOOM
		}
		catch(...){

			#ifdef DEBUG_BLOOM
			cerr << "[" << mpi_rank << "] unable to find valid Bloom parameters" << endl;
			#endif // DEBUG_BLOOM

			// We were unable to find Bloom filter parameters that satisfied the 
			//requested false_positive_probability.

			if(bcount != NULL){

				delete [] bcount;
				bcount = NULL;
			}

			return STATUS_BLOOM_INVALID;
		}
		
		BloomFilter filter(m_param);
		
		filter.unset_all_bits();
		
		BitVector::BLOCK *dst_ptr = filter.ptr();
		const uint64_t num_dst_block = filter.num_block();

		for(size_t h = 0;h < m_param.num_hash;++h){

			BitVector::BLOCK *src_ptr = valid_bits[h].ptr();
			const uint64_t num_src_block = valid_bits[h].num_block();

			for(uint64_t i = 0;i < num_src_block;i += num_dst_block){
				for(uint64_t j = 0;j < num_dst_block;++j){
					dst_ptr[j] |= src_ptr[i + j];
				}
			}
		}
		
		#ifdef DEBUG_BLOOM
		cerr << "[" << mpi_rank << "] Set Bloom filter bits" << endl;
		#endif // DEBUG_BLOOM

		// ** This is the original implementation that is left here as a warning
		// to future generations on how *not* to merge Bloom filters! **
		// The new implementation is better because it:
		//	a) Only accesses two filters at a time (the source and the destination),
		//	   as opposed to the impementation below which touchs *all* of the 
		//	   m_param.num_hash source filters in addition to the destination filter
		//	b) Sets 8-bits at a time by combining BitVector::BLOCK (i.e. byte) elements
		//	   with the binary OR operator. The original implmentation set each bit
		// 	   individually!
		//
		// Set the Bloom filter bits from the counting bloom filter
		//for(size_t i = 0;i < num_count_bloom;++i){
		//
		//	for(size_t h = 0;h < m_param.num_hash;++h){
		//
				// If any of the valid bits (for any of the *valid*
				// hash functions) are set, then set that bit in the
				// output filter.
		//		if( valid_bits[h].get_bit(i) ){
		//
		//			filter.set_bit(i%filter_len);
		//			break;
		//		}
		//	}
		//}
		
		// Clean up the counting filter
		if(bcount != NULL){
			
			delete [] bcount;
			bcount = NULL;
		}

		// Compute the checksum value to safeguard the Bloom filter data only (metadata and Bloom filter
		// parameters are not included in the crc32).
		filter.update_crc32();
					
		filter.set_info(m_info); // Set the metadata
		
		#ifdef DEBUG_BLOOM
			
		cerr << "Bloom filter checksum = " << std::hex <<  filter.get_crc32() << std::dec << endl;
		cerr << "Bloom filter occupancy = " << ( (float)filter.count() )/m_param.filter_len() << endl;

		cerr << "Meta data:" << endl;
		cerr << "\trun_accession: " << accession_to_str(m_info.run_accession) << endl;
		cerr << "\texperiment_accession: " << accession_to_str(m_info.experiment_accession) << endl;
		cerr << "\texperiment_title: " << m_info.experiment_title << endl;
		cerr << "\texperiment_design_description: " << m_info.experiment_design_description << endl;
		cerr << "\texperiment_library_name: " << m_info.experiment_library_name << endl;
		cerr << "\texperiment_library_strategy: " << m_info.experiment_library_strategy << endl;
		cerr << "\texperiment_library_source: " << m_info.experiment_library_source << endl;
		cerr << "\texperiment_library_selection: " << m_info.experiment_library_selection << endl;
		cerr << "\texperiment_instrument_model: " << m_info.experiment_instrument_model << endl;
		cerr << "\tsample_accession: " << accession_to_str(m_info.sample_accession) << endl;
		cerr << "\tsample_taxa: " << m_info.sample_taxa << endl;
		
		if( !m_info.sample_attributes.empty() ){
		
			cerr << "\tsample attributes:" << endl;
			
			for(MAP<string, string>::const_iterator j = m_info.sample_attributes.begin();
				j != m_info.sample_attributes.end();++j){
				
				cerr << "\t\t" << j->first << ": " << j->second << endl;
			}
		}
		
		cerr << "\tstudy_accession: " << accession_to_str(m_info.study_accession) << endl;
		cerr << "\tstudy_title: " << m_info.study_title << endl;
		cerr << "\tstudy_abstract: " << m_info.study_abstract << endl;
		#endif // DEBUG_BLOOM
		
		const string output_file = m_bloom_dir + PATH_SEPARATOR + accession_to_str(m_acc) + ".bloom";

		#ifdef DEBUG_BLOOM
		cerr << "[" << mpi_rank << "] Writing Bloom filter to: " << output_file << endl;
		#endif // DEBUG_BLOOM

		ofstream fout(output_file.c_str(), ios::binary);
		
		if(!fout){
			throw __FILE__ ":main: Unable to open Bloom filter file for writing";
		}
		
		binary_write(fout, filter);
		
		// Make sure that this file is closed and commited to disk
		fout.close();

		#ifdef DEBUG_BLOOM
		cerr << "[" << mpi_rank << "] Finished with Bloom filter construction" << endl;
		#endif // DEBUG_BLOOM
	}
	catch(const char* error){

		if(bcount != NULL){

			delete [] bcount;
			bcount = NULL;
		}

		// Save the error for debugging
		m_progress.error = error;

		#ifdef DEBUG_BLOOM
        cerr << "[" << mpi_rank << "] Caught the error: " << error << endl;
		#endif // DEBUG_BLOOM

        return STATUS_BLOOM_FAIL;
	}
	catch(const std::exception &error){

		if(bcount != NULL){

			delete [] bcount;
			bcount = NULL;
		}

		// Save the error for debugging
		m_progress.error = error.what();

		#ifdef DEBUG_BLOOM
        cerr << "[" << mpi_rank << "] Caught the error: " << error.what() << endl;
		#endif // DEBUG_BLOOM

        return STATUS_BLOOM_FAIL;
	}
	catch(...){

		if(bcount != NULL){
			
			delete [] bcount;
			bcount = NULL;
		}

		#ifdef DEBUG_BLOOM
		cerr << "[" << mpi_rank << "] Caught an unhandled error" << endl;
		#endif // DEBUG_BLOOM

		return STATUS_BLOOM_FAIL;
	}

	return STATUS_BLOOM_SUCCESS;
}

void count_words(CountingBloom *m_count_ptr, vector<BitVector> &m_valid_bits, 
	size_t &m_num_valid_kmer, const ngs::StringRef &m_seq, 
	const size_t &m_hash_seq_mask, const size_t &m_hash_count_mask, 
	const MaestroOptions &m_opt)
{
	// Since the hash_index is used by both the counting Bloom filter and
	// the sequence Bloom filter, it needs to have room to store the 
	// largest number of hash values for each word
	vector<size_t> hash_index( max(NUM_COUNT_HASH, MAX_NUM_HASH) );
	
	#ifdef DEBUG_BLOOM
	const size_t num_count_bloom = m_hash_count_mask + 1;
	#endif // DEBUG_BLOOM

	const char* begin = m_seq.data();
	const char* end = begin + m_seq.size();

	ForEachDuplexWord(begin, end, m_opt.kmer_len)

		if(ValidWord){

			// The counting Bloom filter may have a different length
			// then the sequence Bloom filters, we need to use modulo division
			// to correctly clamp the hash values for the different Bloom
			// filter types.
			bigsi_hash(hash_index, CanonicalWord, m_opt.kmer_len, m_opt.hash_func);

			#ifdef DEBUG_BLOOM
			if( ( (hash_index[0] & m_hash_count_mask) >= num_count_bloom ) ||
				( (hash_index[1] & m_hash_count_mask) >= num_count_bloom ) ||
				( (hash_index[2] & m_hash_count_mask) >= num_count_bloom ) ||
				( (hash_index[3] & m_hash_count_mask) >= num_count_bloom ) ){
				
				cerr << "[" << mpi_rank << "] **bcount overflow!!**" << endl;
				throw __FILE__ ":make_bloom_filter: Buffer overflow!";
			}
			#endif // DEBUG_BLOOM

			// Using a power-of-2 length allows a fast modulo division
			// operation using a bit-mask
			const unsigned char count_first_0 = m_count_ptr[ hash_index[0] & m_hash_count_mask ].first;
			const unsigned char count_first_1 = m_count_ptr[ hash_index[1] & m_hash_count_mask ].first;

			const unsigned char count_second_0 = m_count_ptr[ hash_index[2] & m_hash_count_mask ].second;
			const unsigned char count_second_1 = m_count_ptr[ hash_index[3] & m_hash_count_mask ].second;

			const unsigned char min_count = 
				min( count_first_0, 
					min(count_first_1,
						min(count_second_0, count_second_1) ) );

			// The count is clamped to be <= MAX_COUNT, so don't increment
			// the number of valid kmers or the per-element count once we've
			// reached MAX_COUNT
			if(min_count < m_opt.min_kmer_count){

				if( min_count == (m_opt.min_kmer_count - 1) ){

					++m_num_valid_kmer;

					for(size_t h = 0;h < MAX_NUM_HASH;++h){

						#ifdef DEBUG_BLOOM
						if( (hash_index[h] & m_hash_seq_mask) >= (1ULL << m_opt.max_log_2_filter_len) ){
							
							cerr << "[" << mpi_rank << "] **valid_bits filter overflow!!**" << endl;
							throw __FILE__ ":make_bloom_filter: Bloom filter overflow!";
						}
						#endif // DEBUG_BLOOM
						
						// Using a power-of-2 length allows a fast modulo division
						// operation using a bit-mask
						m_valid_bits[h].set_bit(hash_index[h] & m_hash_seq_mask);
					}
				}

				// Only increment the counting Bloom filter elements 
				// that have a count equal to the min_count.

				// Counting Bloom filter #1 -- two hash functions
				if(count_first_0 == min_count){
					++m_count_ptr[ hash_index[0] & m_hash_count_mask ].first;
				}

				if(count_first_1 == min_count){
					++m_count_ptr[ hash_index[1] & m_hash_count_mask ].first;
				}

				// Counting Bloom filter #2 -- two hash functions
				if(count_second_0 == min_count){
					++m_count_ptr[ hash_index[2] & m_hash_count_mask ].second;
				}

				if(count_second_1 == min_count){
					++m_count_ptr[ hash_index[3] & m_hash_count_mask ].second;
				}

				#ifdef ORIGINAL
				m_count_ptr[ hash_index[0] & m_hash_count_mask ].first = 
					count_first_0 + ( (count_first_0 == min_count) ? 1 : 0 );

				m_count_ptr[ hash_index[1] & m_hash_count_mask ].first = 
					count_first_1 + ( (count_first_1 == min_count) ? 1 : 0 );

				// Counting Bloom filter #2 -- two hash functions
				m_count_ptr[ hash_index[2] & m_hash_count_mask ].second = 
					count_second_0 + ( (count_second_0 == min_count) ? 1 : 0 );
				
				m_count_ptr[ hash_index[3] & m_hash_count_mask ].second = 
					count_second_1 + ( (count_second_1 == min_count) ? 1 : 0 );
				#endif // ORIGINAL
			}
		}

	EndWord
}