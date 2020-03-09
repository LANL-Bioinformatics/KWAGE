#include "ifind.h"

using namespace std;

// Case-insensitive string matching. Return string::npos if the query is not found in
// the subject.
string::size_type ifind(const string &m_subject, const string &m_query)
{
	for(string::const_iterator s = m_subject.begin();s != m_subject.end();++s){
		
		string::const_iterator q = m_query.begin();
		
		for(string::const_iterator i = s;
			( i != m_subject.end() ) && ( q != m_query.end() ) &&
			( tolower(*i) == tolower(*q) );++q,++i){
		}
		
		if( q == m_query.end() ){
			return s - m_subject.begin();
		}
	}
	
	return string::npos;
}
