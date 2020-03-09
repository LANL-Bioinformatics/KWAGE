#ifndef __OUTPUT
#define __OUTPUT

#include <ostream>
#include <string>
#include <iomanip>
#include "bloom.h"

struct MatchResult
{
	unsigned int num_kmers_found;
	unsigned int num_query_kmer;
	FilterInfo subject_info;
	
	MatchResult()
	{
	};
	
	MatchResult(const unsigned int &m_num_kmers_found, 
		const unsigned int &m_num_query_kmer, const FilterInfo &m_subject_info) :
		num_kmers_found(m_num_kmers_found), 
		num_query_kmer(m_num_query_kmer),
		subject_info(m_subject_info)
	{
	};
	
	// The output is sorted in *descending order* by num_kmers_found 
	inline bool operator<(const MatchResult &m_rhs) const
	{
		// Sort in *descending* order by the number of kmer matches
		return (num_kmers_found > m_rhs.num_kmers_found);
	};
};

void write_csv_header(std::ostream &m_out)
{
	m_out << "query,num_kmers,num_kmers_found,percent_kmers_found,sample_metadata\n";
}

template<typename _ITERATOR>
void write_csv(std::ostream &m_out, const std::string &m_query, 
	_ITERATOR m_begin, _ITERATOR m_end)
{
	for(_ITERATOR i = m_begin;i != m_end;++i){
		
		const float norm = i->num_query_kmer ? 1.0f/i->num_query_kmer : 0.0f;
		
		m_out << '"' << m_query << "\"," 
			<< i->num_query_kmer << ',' 
			<< i->num_kmers_found << ','
			<< (100.0f*i->num_kmers_found)*norm << ",\""
			<< i->subject_info.csv_string() << '"' << std::endl;
	}
}

void write_csv_footer(std::ostream &m_out)
{
	// Nothing to do!
}

void write_json_header(std::ostream &m_out, bool m_multiple_query_matches)
{
	if(m_multiple_query_matches){
		m_out << '[';
	}
}

template<typename _ITERATOR>
void write_json(std::ostream &m_out, const std::string &m_query,
	bool m_multiple_query_matches, bool m_first_match,
	const float &m_threshold, _ITERATOR m_begin, _ITERATOR m_end)
{
	const std::string prefix = m_multiple_query_matches ? "\t" : "";
		
	// Format floating point values to have one significant figure after the decimal point
	// and to not use scientific notation (i.e. use fixed notation instead). This is to match
	// the output format of the original bigsi python implementation.
	m_out << ( (m_multiple_query_matches && !m_first_match) ? "," : "") << '\n' << prefix 
		<< "{\n" << prefix << "\t\"query\": \"" << m_query << "\",\n" << prefix 
		<< "\t\"threshold\": " 
		<< std::showpoint << std::setprecision(1) << std::fixed << m_threshold 
		<< ",\n" << prefix << "\t\"results\": [";
	
	for(_ITERATOR i = m_begin;i != m_end;++i){
	
		const float norm = i->num_query_kmer ? 1.0f/i->num_query_kmer : 0.0f;
		
		m_out << ( (i != m_begin) ? "," : "")
			<< "\n" << prefix << "\t\t{\n" << prefix 
				<< "\t\t\t\"percent_kmers_found\": " 
			<< (100.0*i->num_kmers_found)*norm
			<< ",\n" << prefix << "\t\t\t\"num_kmers\": " << i->num_query_kmer 
			<< ",\n" << prefix << "\t\t\t\"num_kmers_found\": " << i->num_kmers_found
			<< ",\n" << prefix << "\t\t\t\"sample_metadata\": {\n" 
			<< i->subject_info.json_string(prefix + "\t\t\t\t")
			<< "\n" << prefix
			<< "\t\t\t}\n" << prefix << "\t\t}";
	}
	
	if(m_begin != m_end){
		m_out << "\n" << prefix << '\t';
	}
	
	m_out << "]\n" << prefix << "}";
}

void write_json_footer(std::ostream &m_out, bool m_multiple_query_matches)
{
	if(m_multiple_query_matches){
		m_out << "\n]\n";
	}
}

#endif // __OUTPUT
