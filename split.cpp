#include "split.h"

using namespace std;

deque<string> split(const string &m_str, char m_delim)
{
	deque<string> ret;

	if( m_str.empty() ){
		return ret;
	}

	ret.push_back( string() );

	for(string::const_iterator i = m_str.begin();i != m_str.end();++i){

		// A '\n' or DOS/Windows '\r' symbol forces the end of the line
		if( (*i == '\r') || (*i == '\n') ){
				break;
		}

		if(*i == m_delim){
				ret.push_back( string() );
		}
		else{
				ret.back().push_back(*i);
		}
	}

	return ret;
}