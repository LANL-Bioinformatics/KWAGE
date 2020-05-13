#ifndef __IFIND
#define __IFIND

#include <string>

// Case-insensitive string matching. Return string::npos if the query is not found in
// the subject.
std::string::size_type ifind(const std::string &m_subject, const std::string &m_query);

#endif // __IFIND
