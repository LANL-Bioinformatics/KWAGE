#include <sstream>
#include <math.h>
#include <zlib.h> // For crc32

#include "bloom.h"

using namespace std;

BloomParam optimal_bloom_param(const uint32_t &m_kmer_len, const size_t &m_num_kmer, const float &m_p,
	const HashFunction &m_func)
{
	// There are some SRA datasets that contain reads that are *less* than the size of a kmer.
	// These datasets will not yeild any valid kmers
	if(m_num_kmer == 0){
		throw __FILE__ ":optimal_bloom_param: No kmers found";
	}
	
	BloomParam ret;
	
	// Even though the optimal bloom filter parameters do not depend on the hash function
	// (we assume that all hash functions are "ideal"), we set the hash function to make sure
	// it gets included in the BloomParam output.
	ret.hash_func = m_func;
	ret.kmer_len = m_kmer_len;

	// These filter parameters are not yet valid
	bool valid = false;
	
	// Perform a grid search to identify the *smallest* Bloom filter 
	// that is less than the specified false positive rate.
	//
	// For Bloom filters of the same length, find the optimal number of hash
	// functions (within the specified search limits).
	for(ret.log_2_filter_len = MIN_LOG_2_BLOOM_FILTER_LEN;
			ret.log_2_filter_len <= MAX_LOG_2_BLOOM_FILTER_LEN;++ret.log_2_filter_len){
		
		float best_p = 10.0f; // Any value > 1.0 to initialize

		// Find the optimal number of hash functions for this filter length
		for(uint32_t num_hash = MIN_NUM_HASH;num_hash <= MAX_NUM_HASH;++num_hash){
			
			const uint32_t len = 1U << ret.log_2_filter_len;
			
			// The per-filter, per-k-mer probability of a false positive
			const double p = pow(1.0 - pow(1.0 - 1.0/len, m_num_kmer*num_hash), num_hash);
			
			if( (p <= m_p) && (p < best_p) ){
				
				best_p = p;
				ret.num_hash = num_hash;
				valid = true;
			}			
		}

		if(valid == true){

			// Since we are searching in order of ascending filter
			// length, return as soon as we satisfy the false 
			// positive requirement.
			return ret;
		}
	}
	
	throw __FILE__ ":optimal_bloom_param: Unable to satisfy Bloom filter probability bound";
	return ret;
}

// Output a limited subset of metadata in CSV
string FilterInfo::csv_string() const
{
	return run_accession;
}

string FilterInfo::json_string(const string &m_prefix) const
{
	stringstream ssout;

	bool wrote_value = false;

	if( !run_accession.empty() ){

		ssout << m_prefix << "\"run\": \"" << run_accession << '"';

		wrote_value = true;
	}

	if( !experiment_accession.empty() ){

		if(wrote_value){
			ssout << ",\n";
		}

		ssout << m_prefix << "\"experiment\": \"" << experiment_accession << '"';

		wrote_value = true;
	}

	if( !experiment_title.empty() ){

		if(wrote_value){
			ssout << ",\n";
		}

		ssout << m_prefix << "\"experiment title\": \"" << experiment_title << '"';

		wrote_value = true;
	}

	if( !experiment_title.empty() ){

		if(wrote_value){
			ssout << ",\n";
		}

		ssout << m_prefix << "\"experiment title\": \"" << experiment_title << '"';

		wrote_value = true;
	}

	if( !experiment_design_description.empty() ){

		if(wrote_value){
			ssout << ",\n";
		}

		ssout << m_prefix << "\"experiment design\": \"" << experiment_design_description << '"';

		wrote_value = true;
	}

	if( !experiment_library_name.empty() ){

		if(wrote_value){
			ssout << ",\n";
		}

		ssout << m_prefix << "\"experiment library name\": \"" << experiment_library_name << '"';

		wrote_value = true;
	}

	if( !experiment_library_strategy.empty() ){

		if(wrote_value){
			ssout << ",\n";
		}

		ssout << m_prefix << "\"experiment library strategy\": \"" << experiment_library_strategy << '"';

		wrote_value = true;
	}

	if( !experiment_library_source.empty() ){

		if(wrote_value){
			ssout << ",\n";
		}

		ssout << m_prefix << "\"experiment library source\": \"" << experiment_library_source << '"';

		wrote_value = true;
	}

	if( !experiment_library_selection.empty() ){

		if(wrote_value){
			ssout << ",\n";
		}

		ssout << m_prefix << "\"experiment library selection\": \"" << experiment_library_selection << '"';

		wrote_value = true;
	}

	if( !experiment_instrument_model.empty() ){

		if(wrote_value){
			ssout << ",\n";
		}

		ssout << m_prefix << "\"experiment instrument model\": \"" << experiment_instrument_model << '"';

		wrote_value = true;
	}

	if( !sample_accession.empty() ){

		if(wrote_value){
			ssout << ",\n";
		}

		ssout << m_prefix << "\"sample\": \"" << sample_accession << '"';

		wrote_value = true;
	}

	if( !sample_taxa.empty() ){

		if(wrote_value){
			ssout << ",\n";
		}

		ssout << m_prefix << "\"sample taxa\": \"" << sample_taxa << '"';

		wrote_value = true;
	}

	if( !sample_attributes.empty() ){

		if(wrote_value){
			ssout << ",\n";
		}

		bool wrote_attribute_value = false;

		ssout << m_prefix << "\"sample attributes\": [\n";

		for(MULTIMAP<string, string>::const_iterator i = sample_attributes.begin();i != sample_attributes.end();++i){

			if(wrote_attribute_value){
				ssout << ",\n";
			}

			ssout << m_prefix << "\t{\n";
			ssout << m_prefix << "\t\t\"tag\": \"" << i->first << "\",\n";
			ssout << m_prefix << "\t\t\"value\": \"" << i->second << "\"\n";
			ssout << m_prefix << "\t}";

			wrote_attribute_value = true;
		}

		ssout <<'\n' << m_prefix << ']';

		wrote_value = true;
	}

	if( !study_accession.empty() ){

		if(wrote_value){
			ssout << ",\n";
		}

		ssout << m_prefix << "\"study\": \"" << study_accession << '"';

		wrote_value = true;
	}

	if( !study_title.empty() ){

		if(wrote_value){
			ssout << ",\n";
		}

		ssout << m_prefix << "\"study title\": \"" << study_title << '"';

		wrote_value = true;
	}

	if( !study_abstract.empty() ){

		if(wrote_value){
			ssout << ",\n";
		}

		ssout << m_prefix << "\"study abstract\": \"" << study_abstract << '"';

		wrote_value = true;
	}

	return ssout.str();
}

unsigned int BitVector::crc32() const
{
	// Use the crc32 function provided by zlib
	return ::crc32( 0U, buffer, num_block() );
};

unsigned int BloomFilter::update_crc32()
{
	bloom_crc32 = BitVector::crc32();
	
	return bloom_crc32;
}

bool BloomFilter::test_crc32() const
{
	return ( bloom_crc32 == BitVector::crc32() );
}
