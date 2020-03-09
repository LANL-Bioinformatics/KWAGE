#ifndef __BINARY_IO
#define	__BINARY_IO

#include <ostream>
#include <istream>
#include <deque>
#include <vector>
#include <string>
#include <unordered_map>

// Forward function definitions for containers (needed to be able to transport nested C++ structures):
template<class T> void binary_write(std::ostream &m_out, const std::deque<T> &m_obj);
template<class T> void binary_read(std::istream &m_in, std::deque<T> &m_obj);

template<class T> void binary_write(std::ostream &m_out, const std::vector<T> &m_obj);
template<class T> void binary_read(std::istream &m_in, std::vector<T> &m_obj);

template<class A, class B> void binary_write(std::ostream &m_out, const std::pair<A, B> &m_obj);
template<class A, class B> void binary_read(std::istream &m_in, std::pair<A, B> &m_obj);

template<class A, class B> void binary_write(std::ostream &m_out, const std::unordered_map<A,B> &m_obj);
template<class A, class B> void binary_read(std::istream &m_in, std::unordered_map<A,B> &m_obj);

template<class A, class B> void binary_write(std::ostream &m_out, const std::unordered_multimap<A, B> &m_obj);
template<class A, class B> void binary_read(std::istream &m_in, std::unordered_multimap<A, B> &m_obj);

// Use a template for simple objects. Specialize as needed for more complex types
template <class T>
void binary_write(std::ostream &m_out, const T &m_obj)
{
	// Force a *compile* time test of whether this is a native or derived type
	static_assert(std::is_fundamental<T>::value || std::is_enum<T>::value,
		":mpi_pack: Non-fundamental or non-enum type passed as template");
	
	m_out.write( (char*)&m_obj, sizeof(m_obj) );
}

template <class T>
void binary_read(std::istream &m_in, T &m_obj)
{
	// Force a *compile* time test of whether this is a native or derived type
	static_assert(std::is_fundamental<T>::value || std::is_enum<T>::value,
		":mpi_unpack: Non-fundamental or non-enum type passed as template");
	
	m_in.read( (char*)&m_obj, sizeof(m_obj) );
}

// Specialization for string
template<>
	void binary_write(std::ostream &m_out, const std::string &m_str);
template<>
	void binary_read(std::istream &m_in, std::string &m_str);

// Specialization for __uint128_t (which is not recognized by g++ as either a fundamental
// integral or scalar type!). __uint128_t is used to store both accessions and taxonomy ids ...
template<>
	void binary_write(std::ostream &m_out, const __uint128_t &m_obj);
template<>
	void binary_read(std::istream &m_in, __uint128_t &m_obj);
	
/////////////////////////////////////////////////////////////////////////////////////////
// Overload for std::deque
/////////////////////////////////////////////////////////////////////////////////////////

template<class T>
void binary_write(std::ostream &m_out, const std::deque<T> &m_obj)
{
	binary_write(m_out, m_obj.size() );

	for(typename std::deque<T>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		binary_write(m_out, *i);
	}
}

template<class T>
void binary_read(std::istream &m_in, std::deque<T> &m_obj)
{
	size_t len;
	
	binary_read(m_in, len);
	
	m_obj.resize(len);
	
	for(size_t i = 0;i < len;++i){
		binary_read(m_in, m_obj[i]);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// Overload for std::vector
/////////////////////////////////////////////////////////////////////////////////////////

template<class T>
void binary_write(std::ostream &m_out, const std::vector<T> &m_obj)
{
	binary_write( m_out, m_obj.size() );
	
	for(typename std::vector<T>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		binary_write(m_out, *i);
	}
}

template<class T>
void binary_read(std::istream &m_in, std::vector<T> &m_obj)
{
	size_t len;
	
	binary_read(m_in, len);
	
	m_obj.resize(len);
	
	for(size_t i = 0;i < len;++i){
		binary_read(m_in, m_obj[i]);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// Overload for std::pair
/////////////////////////////////////////////////////////////////////////////////////////

template<class A, class B>
void binary_write(std::ostream &m_out, const std::pair<A, B> &m_obj)
{
	binary_write(m_out, m_obj.first);
	binary_write(m_out, m_obj.second);
}

template<class A, class B>
void binary_read(std::istream &m_in, std::pair<A, B> &m_obj)
{       
    binary_read(m_in, m_obj.first);
	binary_read(m_in, m_obj.second);
}

/////////////////////////////////////////////////////////////////////////////////////////
// Overload for std::unordered_map
/////////////////////////////////////////////////////////////////////////////////////////

template<class A, class B>
void binary_write(std::ostream &m_out, const std::unordered_map<A,B> &m_obj)
{
	binary_write( m_out, m_obj.size() );
	
	for(typename std::unordered_map<A,B>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		
		binary_write(m_out, i->first);
		binary_write(m_out, i->second);
	}
}

template<class A, class B>
void binary_read(std::istream &m_in, std::unordered_map<A,B> &m_obj)
{
	size_t len;
	
	binary_read(m_in, len);

	m_obj.clear();
	
	for(size_t i = 0;i < len;++i){
		
		std::pair<A,B> local;
		
		binary_read(m_in, local.first);
		binary_read(m_in, local.second);
		
		m_obj.insert(local);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// Overload for std::unordered_multimap
/////////////////////////////////////////////////////////////////////////////////////////

template<class A, class B>
void binary_write(std::ostream &m_out, const std::unordered_multimap<A, B> &m_obj)
{	
	binary_write( m_out, m_obj.size() );

	for(typename std::unordered_multimap<A, B>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		
		binary_write(m_out, i->first);
		binary_write(m_out, i->second);
	}
}

template<class A, class B>
void binary_read(std::istream &m_in, std::unordered_multimap<A, B> &m_obj)
{
	size_t len;
	
	binary_read(m_in, len);
	
	m_obj.clear();
	
	for(size_t i = 0;i < len;++i){
		
		std::pair<A, B> local;
		
		binary_read(m_in, local.first);
		binary_read(m_in, local.second);
		
		m_obj.insert(local);
	}
}

#endif // __BINARY_IO
