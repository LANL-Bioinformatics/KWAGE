// Search a KWAGE Bloom filter database using input fasta sequences as queries
// J. D. Gans
// Bioscience Division, B-10
// Los Alamos National Laboratory
// Mon Oct 28 17:42:59 2019

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdlib.h>

#include "kwage.h"
#include "options.h"
#include "word.h"
#include "sort.h"
#include "bloom.h"
#include "hash.h"
#include "output.h"
#include "parse_sequence.h"
#include "keys.h"

using namespace std;

// Return true if the specified query matched one or more subject Bloom filters
bool search(unordered_map< size_t, deque<MatchResult> > &m_search_results, 
	ifstream &m_fsubject, const string &m_query, const size_t &m_query_id,
	const unsigned long int &m_bloom_start, const unsigned long int &m_info_start,
	const size_t &m_slice_size,
	const DBFileHeader &m_header,
	const SearchOptions &m_opt);

// Global variables for MPI
int mpi_numtasks;
int mpi_rank;

int main(int argc, char *argv[])
{
	try{
		time_t profile = time(NULL);
		
		SearchOptions opt(argc, argv);
		
		if(opt.quit){
			return EXIT_SUCCESS;
		}
		
		ofstream fout;
		
		if( !opt.output_file.empty() ){
		
			fout.open( opt.output_file.c_str() );
			
			if(!fout){
				
				cerr << "Unable to open " << opt.output_file << " for writing" << endl;
				return EXIT_FAILURE;
			}
		}
		
		ostream &out = fout.is_open() ? fout : cout;
		
		// Store all of the results in memory (to allow rank 0 to collect and sort
		// all of the search results for the user).
		unordered_map< size_t /*query id*/, deque<MatchResult> > file_search_results;
		unordered_map< size_t /*query id*/, deque<MatchResult> > command_line_search_results;
		
		// We will store the file-based query deflines that match one or more subject Bloom filters
		unordered_map<size_t /*query id*/, string /*defline*/> file_query_info;
		
		const size_t num_subject_files = opt.subject_files.size();
		
		// Multithreading does not offer much (if any) advantage when searching uncompressed
		// database files. However, it does effectively hide the cost of inflating compressed
		// bit slice buffers.
		#pragma omp parallel
		{
			try{

			// Local copies of the search results to allow individual threads to
			// to progress asyncronously
			unordered_map< size_t, deque<MatchResult> > local_file_search_results;
			unordered_map< size_t, deque<MatchResult> > local_command_line_search_results;
			unordered_map<size_t, string /*defline*/> local_file_query_info;
		
			#pragma omp for
			for(size_t file_index = 0;file_index < num_subject_files;++file_index){

				ifstream fsubject(opt.subject_files[file_index].c_str(), ios::binary);

				if(!fsubject){

					cerr << "Unable to open database file " << opt.subject_files[file_index]
						<< " for reading" << endl;
					throw __FILE__ ":main: I/O error";
				}

				// Read the database header
				DBFileHeader header;

				binary_read(fsubject, header);

				if(!fsubject){
					throw __FILE__ ":main: Unable to read header";
				}

				// The size of each uncompressed bitslice in bytes
				const size_t slice_size = header.num_filter/BloomFilter::BITS_PER_BLOCK + 
					( (header.num_filter%BloomFilter::BITS_PER_BLOCK == 0) ? 0 : 1 );

				// Mark the start of the bit-sliced bloom filter data
				const unsigned long int bloom_start = fsubject.tellg();
				const unsigned long int info_start = header.info_start;

				// Search sequences provided on the command line
				for(deque<string>::const_iterator query_iter = opt.query_seq.begin();
					query_iter != opt.query_seq.end();++query_iter){

					search(local_command_line_search_results, fsubject, *query_iter, 
						query_iter - opt.query_seq.begin(), // The query id
						bloom_start, info_start,
						slice_size, header,
						opt);
				}

				// Search sequences provided in a sequence file
				size_t query_id = 0;

				for(deque<string>::const_iterator query_file_iter = opt.query_files.begin();
					query_file_iter != opt.query_files.end();++query_file_iter){

					// Iterate through the contents of a sequence file (currently fasta or fastq)
					SequenceIterator seq_iter(*query_file_iter);

					while(seq_iter){
						
						if( search(local_file_search_results, fsubject, seq_iter.get_seq(), query_id,
							bloom_start, info_start,
							slice_size, header,
							opt) ){

							local_file_query_info[query_id] = seq_iter.get_info();
						}

						++seq_iter;
						++query_id;
					}
				}

				fsubject.close();
			}
			
			// Merge the local (i.e. per-thread) results into the global results
			#pragma omp critical
			{
				for(unordered_map< size_t, deque<MatchResult> >::const_iterator i = local_file_search_results.begin();
					i != local_file_search_results.end();++i){

					deque<MatchResult> &ref = file_search_results[i->first];

					ref.insert( ref.end(), i->second.begin(), i->second.end() );
				}

				for(unordered_map< size_t, deque<MatchResult> >::const_iterator i = local_command_line_search_results.begin();
					i != local_command_line_search_results.end();++i){

					deque<MatchResult> &ref = command_line_search_results[i->first];

					ref.insert( ref.end(), i->second.begin(), i->second.end() );
				}

				for(unordered_map<size_t, string>::const_iterator i = local_file_query_info.begin();
					i != local_file_query_info.end();++i){
					
					file_query_info[i->first] = i->second;
				}
			}

			}
			catch(const char *error){
				cerr << "Caught the search error: " << error << endl;
				throw error;
			}
			catch(...){
				cerr << "Caught an unhandled search error" << endl;
				throw "Unhandled search error";
			}
		}
		
		// Sort all of the results
		for(unordered_map< size_t, deque<MatchResult> >::iterator i = command_line_search_results.begin();
			i != command_line_search_results.end();++i){
			
			SORT( i->second.begin(), i->second.end() );
		}
		
		for(unordered_map< size_t, deque<MatchResult> >::iterator i = file_search_results.begin();
			i != file_search_results.end();++i){
			
			SORT( i->second.begin(), i->second.end() );
		}
		
		// Did we match multiple queries? If so, we will need to add a header, footer
		// and indentation to any JSON file.
		const bool multiple_query_matches = 
			( command_line_search_results.size() + file_search_results.size() ) > 1;
		
		switch(opt.output_format){
			case SearchOptions::OUTPUT_CSV:
			
				write_csv_header(out);
				break;
			case SearchOptions::OUTPUT_JSON:
				
				write_json_header(out, multiple_query_matches);
				break;
			default:
				throw __FILE__ ":main: Unknown output file format (1)";
		};

		bool first_match = true;
				
		// Print any command line matches first
		vector<size_t> id = keys(command_line_search_results);
		
		SORT( id.begin(), id.end() );

		for(vector<size_t>::const_iterator i = id.begin();i != id.end();++i){
			
			unordered_map< size_t, deque<MatchResult> >::const_iterator iter = 
				command_line_search_results.find(*i);
				
			if( iter == command_line_search_results.end() ){
				throw __FILE__ ":main: Unable to lookup query id in command_line_search_results";
			}
			
			stringstream ssin;
			
			// Make dummy information for command line sequences
			ssin << "command line seq " << *i;
			
			switch(opt.output_format){
				case SearchOptions::OUTPUT_CSV:

					write_csv( out, ssin.str(), 
						iter->second.begin(), iter->second.end() );

					break;
				case SearchOptions::OUTPUT_JSON:

					write_json( out, ssin.str(),
						multiple_query_matches, first_match,
						opt.threshold,
						iter->second.begin(), iter->second.end() );
					break;
				default:
					throw __FILE__ ":main: Unknown output file format (2)";
			};
			
			first_match = false;
		}
		
		// Print the file-based matches
		id = keys(file_search_results);
		
		SORT( id.begin(), id.end() );

		for(vector<size_t>::const_iterator i = id.begin();i != id.end();++i){
			
			unordered_map< size_t, deque<MatchResult> >::const_iterator iter = 
				file_search_results.find(*i);
				
			if( iter == file_search_results.end() ){
				throw __FILE__ ":main: Unable to lookup query id in file_search_results";
			}
			
			unordered_map<size_t, string>::const_iterator info_iter = file_query_info.find(*i);
			
			if( info_iter == file_query_info.end() ){
				throw __FILE__ ":main: Unable to lookup query id in file_query_info";
			}
			
			switch(opt.output_format){
				case SearchOptions::OUTPUT_CSV:

					write_csv( out, info_iter->second, 
						iter->second.begin(), iter->second.end() );

					break;
				case SearchOptions::OUTPUT_JSON:

					write_json( out, info_iter->second,
						multiple_query_matches, first_match,
						opt.threshold,
						iter->second.begin(), iter->second.end() );
					break;
				default:
					throw __FILE__ ":main: Unknown output file format (2)";
			};
			
			first_match = false;
		}
		
		switch(opt.output_format){
			case SearchOptions::OUTPUT_CSV:
			
				write_csv_footer(out);
				break;
			case SearchOptions::OUTPUT_JSON:
				
				write_json_footer(out, multiple_query_matches);
				break;
			default:
				throw __FILE__ ":main: Unknown output file format (1)";
		};
		
		profile = time(NULL) - profile;
		
		cerr << "Search complete in " << profile << " sec" << endl;
		
	}
	catch(const char *error){
		cerr << "Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(const string error){
		cerr << "Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(...){
		cerr << "Caught an unhandled error" << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}

// Return true if we found one or more matches for the specifed query,
// false otherwise
bool search(unordered_map< size_t, deque<MatchResult> > &m_search_results, 
	ifstream &m_fsubject, const string &m_query, const size_t &m_query_id,
	const unsigned long int &m_bloom_start, const unsigned long int &m_info_start,
	const size_t &m_slice_size,
	const DBFileHeader &m_header,
	const SearchOptions &m_opt)
{
	bool ret = false;
	
	const bool complete_match = (m_opt.threshold == 1.0f);
	
	// Version 1: Extract the set of query kmers
	deque<Word> kmers;

	ForEachDuplexWord(m_query.c_str(), m_query.c_str() + m_query.size(), m_header.kmer_len)

		if(ValidWord){
			kmers.push_back(CanonicalWord);
		}

	EndWord

	SORT( kmers.begin(), kmers.end() );

	deque<Word>::const_iterator end_kmers = unique( kmers.begin(), kmers.end() );

	const unsigned int num_query_kmer = end_kmers - kmers.begin();
	
	// The query sequence is too small for the current kmer size
	if(num_query_kmer == 0){
		return false;
	}

	BitVector complete_match_mask;
	vector<unsigned int> match_count;
	unsigned int query_threshold = 0;

	if(complete_match){

		complete_match_mask.resize(m_header.num_filter);

		// Initialize the match bit array to "true" for every Bloom filter. 
		// We will reduce this number by successive binary AND operations
		complete_match_mask.set_all_bits();
	}
	else{					
		match_count.resize(m_header.num_filter);

		query_threshold = m_opt.threshold*num_query_kmer;
	}
	
	BitVector kmer_match(m_header.num_filter);
	BitVector slice(m_header.num_filter);
	
	// The worst case search scenario (in terms of terminating the search early) is
	// no match to any Bloom filter from kmer 0 to kmer (1.0f - m_opt.threshold)*num_query_kmer and
	// then matches for the remainder of the kmers
	deque<Word>::const_iterator mid_kmers = kmers.begin() + (1.0f - m_opt.threshold)*num_query_kmer;
	
	const size_t filter_len = m_header.filter_len();
	
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Test for matches to the first (1.0f - m_opt.threshold)*num_query_kmer kmers *without* attempting
	// to exit early.
	for(deque<Word>::const_iterator i = kmers.begin();i != mid_kmers;++i){

		kmer_match.set_all_bits();

		// The BIGSI python implementation starts from 0 when seeding the hash algorithm
		for(size_t h = 0;h < m_header.num_hash;++h){

			const size_t slice_index = bigsi_hash(*i, m_header.kmer_len, h, 
				m_header.hash_func)%filter_len;
			
			m_fsubject.seekg(m_bloom_start + slice_index*m_slice_size);
				
			slice.read(m_fsubject);
				
			if(!m_fsubject){
				throw __FILE__":search: Error reading slice from file (1)";
			}

			kmer_match &= slice;
		}

		// Separate code paths for threshold == 1.0 (bitmask) 
		// and threshold < 1.0 (vector)
		if(complete_match){
			complete_match_mask &= kmer_match;
		}
		else{
			kmer_match.increment_count(match_count);
		}
	}
	
	/////////////////////////////////////////////////////////////////////////////////
	// Test for matches to the remaining kmers and check to see if we can exit early
	for(deque<Word>::const_iterator i = mid_kmers;i != end_kmers;++i){

		kmer_match.set_all_bits();

		// The BIGSI python implementation starts from 0 when seeding the hash algorithm
		for(size_t h = 0;h < m_header.num_hash;++h){

			const size_t slice_index = bigsi_hash(*i, m_header.kmer_len, h, 
				m_header.hash_func)%filter_len;
			
			m_fsubject.seekg(m_bloom_start + slice_index*m_slice_size);

			slice.read(m_fsubject);
			
			if(!m_fsubject){
				throw __FILE__":search: Error reading slice from file (2)";
			}

			kmer_match &= slice;
		}

		// Separate code paths for threshold == 1.0 (bitmask) 
		// and threshold < 1.0 (vector)
		if(complete_match){
			
			complete_match_mask &= kmer_match;
			
			// Breakout early if we cannot find at least one 
			// Bloom filter that matches all kmers
			if( !complete_match_mask.max_bit() ){
				
				// If no bit is set, we can stop searching
				break;
			}
		}
		else{
			
			kmer_match.increment_count(match_count);
			
			// Breakout early if even the best matching Bloom filter 
			// does not have enough matches
			if( ( *max_element( match_count.begin(), match_count.end() ) + 
				(end_kmers - i) ) < query_threshold){
				break;
			}
		}
	}
	
	// Store the match results for the current query
	unordered_map< size_t, deque<MatchResult> >::iterator result_iter = m_search_results.end();
	
	// Is this query a match to one or more target sequences?
	for(size_t i = 0;i < m_header.num_filter;++i){

		bool matched_filter = false;

		if(complete_match){
			matched_filter = complete_match_mask.get_bit(i);
		}
		else{
			matched_filter = (match_count[i] >= query_threshold);
		}

		if(matched_filter){

			// Read the Bloom filter info id in two steps:
			// 1) Read the address of the group id
			// 2) Go to this address to read the group id string
			m_fsubject.seekg( m_info_start + i*sizeof(unsigned long int) );

			unsigned long int info_loc;

			m_fsubject.read( (char*)&info_loc, sizeof(unsigned long int) );

			m_fsubject.seekg(info_loc);

			FilterInfo info;

			binary_read(m_fsubject, info);

			const unsigned int num_match = complete_match ? 
				num_query_kmer : match_count[i];
			
			// Is this the first match for this query against any Bloom filter in the
			// current file?
			if( result_iter == m_search_results.end() ){
				
				result_iter = m_search_results.find(m_query_id);
				
				// If this is the first match to a query for *any* file, then we will
				// need to initialize the search results to make room for the match record.
				if( result_iter == m_search_results.end() ){
					result_iter = m_search_results.insert( 
						make_pair(m_query_id, deque<MatchResult>() ) ).first;
				}
			}
			
			result_iter->second.push_back( MatchResult(num_match, num_query_kmer, info) );
			
			ret = true;
		}
	}
	
	return ret;
}
