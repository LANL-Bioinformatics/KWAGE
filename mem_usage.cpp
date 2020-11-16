#include "mem_usage.h"
#include "string_conversion.h"

#include <fstream>

using namespace std;

size_t parse_meminfo(const string &m_line);

// Return the % of total memory used
float memory_usage()
{
	ifstream fin("/proc/meminfo");

	if(!fin){
		return false;
	}

	double total = 0.0;
	double avail = 0.0;

	string line;
	bool read_total = false;
	bool read_avail = false;

	while( !(read_total && read_avail) && getline(fin, line) ){

		if(line.find("MemTotal:") != string::npos){

			total = parse_meminfo(line);
			read_total = true;
		}

		if(line.find("MemAvailable:") != string::npos){

			avail = parse_meminfo(line);
			read_avail = true;
		}
	}

	if(total <= 0.0){
		return -1.0f;
	}

	return 100*(total - avail)/total;
}

size_t parse_meminfo(const string &m_line)
{
	// Return the second column of the input line
	string::const_iterator begin = m_line.begin();

	// Skip over the key name
	while( ( begin != m_line.end() ) && !isspace(*begin) ){
		++begin;
	}

	// Skip over the white space
	while( ( begin != m_line.end() ) && isspace(*begin) ){
		++begin;
	}

	// This is value
	string::const_iterator end = begin + 1;

	// Find the end of the value
	while( ( end != m_line.end() ) && !isspace(*end) ){
		++end;
	}

	if(begin >= end){
		throw __FILE__ ":parse_meminfo: Unable to parse value";
	}

	return str_to_size_t( string(begin, end) );
}