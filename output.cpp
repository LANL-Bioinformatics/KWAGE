#include "output.h"
#include <iostream>
#include <iomanip>

using namespace std;

// The origianl bigsi python implementation uses four spaces instead of a tab
#define		TAB	"    "

void write_csv_header(ostream &m_out)
{
	m_out << "\"query\",\"num_kmers\",\"num_kmers_found\",\"percent_kmers_found\",\"sample_name\"" << endl;
}

void write_json_header(ostream &m_out)
{
	// No global header
}

void write_json_query_header(ostream &m_out, const string &m_query, const float &m_threshold)
{
	// Format floating point values to have one significant figure after the decimal point
	// and to not use scientific notation (i.e. use fixed notation instead). This is to match
	// the output format of the original bigsi python implementation.
	m_out << "{\n" << TAB << "\"query\": \"" << m_query << "\",\n" TAB "\"threshold\": " 
		<< std::showpoint << std::setprecision(1) << std::fixed << m_threshold 
		<< ",\n" TAB "\"results\": [";
}


void write_csv_match(ostream &m_out, const string &m_query, const unsigned int &m_num_query_kmer, 
	const unsigned int &m_num_match_kmer, const string &m_group_id)
{
	m_out << '"' << m_query << "\",\"" 
		<< m_num_query_kmer << ',' 
		<< m_num_match_kmer << ','
		<< (100.0f*m_num_match_kmer)/m_num_query_kmer << ",\""
		<< m_group_id << '"' << endl;
}

void write_json_match(ostream &m_out, const unsigned int &m_num_query_kmer, 
	const unsigned int &m_num_match_kmer, bool m_found_match, const string &m_group_id)
{
	m_out << (m_found_match ? "," : "")
		<< "\n" TAB TAB "{\n" TAB TAB TAB "\"percent_kmers_found\": " 
		<< (100.0*m_num_match_kmer)/m_num_query_kmer 
		<< ",\n" TAB TAB TAB "\"num_kmers\": " << m_num_query_kmer 
		<< ",\n" TAB TAB TAB "\"num_kmers_found\": " << m_num_match_kmer
		<< ",\n" TAB TAB TAB "\"sample_name\": \"" << m_group_id 
		<< "\"\n" TAB TAB "}";
}

void write_json_query_footer(ostream &m_out, bool m_found_match)
{
	m_out << (m_found_match ? "\n" TAB : "")
		<< "],\n" TAB "\"citation\": \"http://dx.doi.org/10.1038/s41587-018-0010-1\"\n}" << endl;
}
