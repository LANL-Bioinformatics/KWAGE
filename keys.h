#ifndef __KEYS
#define __KEYS

#include <unordered_map>
#include <vector>
#include <deque>
#include <algorithm>

#include "sort.h"

template<class A, class B>
inline std::vector<A> keys(const std::unordered_multimap<A, B> &m_db)
{
	// Accumulate with a deque
	std::deque<A> ret;

	if( m_db.empty() ){
		return std::vector<A>();
	}
	
	ret.push_back(m_db.begin()->first);
	
	for(typename std::unordered_multimap<A, B>::const_iterator i = m_db.begin();i != m_db.end();++i){
	
		// Can't use x != y, so use !(x==y)
		if( !(ret.back() == i->first) ){
			ret.push_back(i->first);
		}
	}

	// Use the gnu parallel sort 
	SORT( ret.begin(), ret.end() );

	// Return the unique elements as a vector
	return std::vector<A>( ret.begin(),  std::unique( ret.begin(), ret.end() ) );
};

template<class A, class B>
inline std::vector<A> keys(const std::unordered_map<A, B> &m_db)
{
	// Accumulate with a deque
	std::deque<A> ret;

	if( m_db.empty() ){
		return std::vector<A>();
	}
	
	ret.push_back(m_db.begin()->first);
	
	for(typename std::unordered_map<A, B>::const_iterator i = m_db.begin();i != m_db.end();++i){
	
		// Can't use x != y, so use !(x==y)
		if( !(ret.back() == i->first) ){
			ret.push_back(i->first);
		}
	}

	// Use the gnu parallel sort 
	SORT( ret.begin(), ret.end() );

	// Return the unique elements as a vector
	return std::vector<A>( ret.begin(),  std::unique( ret.begin(), ret.end() ) );
};

#endif // __KEYS
