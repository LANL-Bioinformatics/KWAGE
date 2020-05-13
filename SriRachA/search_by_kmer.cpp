#include "sriracha.h"

#include <iostream>

#include "word.h"
#include "sort.h"
#include "options.h"

using namespace std;

extern int mpi_rank;
extern int mpi_numtasks;

void search_by_kmer(const string &m_seq, const unsigned int &m_read_index, const unsigned int &m_read_subindex, void* m_param[])
{
	// Validate the input parameters
	if(m_param == NULL){
		throw __FILE__ ":search_by_kmer: m_param == NULL";
	}

	if(m_param[0] == NULL){
		throw __FILE__ ":search_by_kmer: m_param[0] == NULL (results)";
	}

	vector< deque<SearchMatch> > *results_ptr = (vector< deque<SearchMatch> > *)(m_param[0]);

	if(m_param[1] == NULL){
		throw __FILE__ ":search_by_kmer: m_param[1] == NULL (subject_kmers)";
	}

	deque< pair< string /*defline*/, deque<Word> /*kmers*/> > *subject_kmers_ptr = 
		(deque< pair< string, deque<Word> > >*)(m_param[1]);

	if(m_param[2] == NULL){
		throw __FILE__ ":search_by_kmer: m_param[2] == NULL (num_perfect_match)";
	}

	vector<size_t> *num_perfect_match_ptr = (vector<size_t> *)(m_param[2]);

	if(m_param[3] == NULL){
		throw __FILE__ ":search_by_kmer: m_param[2] == NULL (opt)";
	}

	SrirachaOptions *opt_ptr = (SrirachaOptions *)(m_param[3]);

	if(opt_ptr->verbose >= CHATTY){
		cerr << "[" << mpi_rank << "] is searching read " << m_read_index << endl; 
	}

	if(m_seq.size() < opt_ptr->min_read_length){
		return;
	}

	deque<Word> local_kmers;

	ForEachDuplexWord(m_seq, opt_ptr->kmer_len)

		if(ValidWord){
			local_kmers.push_back(CanonicalWord);
		}

	EndWord

	const size_t num_kmer = local_kmers.size();

	// Count the number of valid bases
	if( num_kmer < opt_ptr->min_valid_kmer ){
		return;
	}

	// Turn the list of extracted kmers into a set
	SORT( local_kmers.begin(), local_kmers.end() );
	local_kmers.erase( unique( local_kmers.begin(), local_kmers.end() ), local_kmers.end() );

	const size_t num_unique_kmer = local_kmers.size();

	// Make sure that this read has sufficient complexity
	if( float(num_unique_kmer)/num_kmer < opt_ptr->min_read_complexity ){
		return;
	}

	const size_t match_threshold = max( size_t(1), size_t(num_unique_kmer * opt_ptr->kmer_match_threshold) );

	const size_t num_subject = subject_kmers_ptr->size();

	for(size_t index = 0;index < num_subject;++index){

		// If we have already found the maximum number of perfect matches,
		// we don't need to keep searching.
		if( (*num_perfect_match_ptr)[index] >= opt_ptr->max_num_match ){
			continue;
		}

		const deque<Word> &subject_kmers = (*subject_kmers_ptr)[index].second;

		// Compute the size of the intersection between the query and subject sets. Note that I directly
		// compared the approach below to the approach of directly iterating through the query and subject
		// sets. On a mid-sized AWS instance, the approach used below is almost 1.6-fold faster. The test set
		// was Illumina datasets (short, ~150 bp reads), and it may be worth testing a long read data
		// set (like PacBio or Oxford Nanopore).
		size_t count = 0;

		for(deque<Word>::const_iterator q = local_kmers.begin();q != local_kmers.end();++q ){

			// Is this a match?
			deque<Word>::const_iterator iter = lower_bound( subject_kmers.begin(), subject_kmers.end(), *q);

			count += ( ( iter != subject_kmers.end() ) && (*iter == *q) ) ? 1 : 0;
		}
	
		if(count >= match_threshold){
			
			deque<SearchMatch> &results = (*results_ptr)[index];

			const float score = float(count)/num_unique_kmer;

			results.push_back( SearchMatch(m_read_index, m_read_subindex, score) );

			if(score == 1.0){
				++( (*num_perfect_match_ptr)[index] );
			}

			if( (opt_ptr->max_num_match > 0) && ( results.size() > 10*opt_ptr->max_num_match ) ){

				// Cull the matches to avoid running out of RAM
				SORT( results.begin(), results.end() );

				results.resize(opt_ptr->max_num_match);
			}
		}
	}
}