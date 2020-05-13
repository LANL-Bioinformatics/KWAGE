#include "mpi_util.h"
#include "options.h"
#include "date.h"
#include "sriracha.h"
#include <string.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for std::string
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<string>(const string &m_str)
{
	return sizeof(size_t) + m_str.size();
}

template<>
unsigned char* mpi_pack<string>(unsigned char* m_ptr, const string &m_str)
{
	size_t len = m_str.size();
	
	memcpy( m_ptr, &len, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	memcpy(m_ptr, m_str.c_str(), len);
	m_ptr += len;
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack<string>(unsigned char* m_ptr, string &m_str)
{
	size_t len;
	
	memcpy( &len, m_ptr, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	m_str.assign( (char*)m_ptr, len );
	m_ptr += len;
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for __uint128_t
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<__uint128_t>(const __uint128_t &m_obj)
{
	return sizeof(__uint128_t);
}

template<>
unsigned char* mpi_pack<__uint128_t>(unsigned char* m_ptr, const __uint128_t &m_obj)
{
	memcpy( m_ptr, &m_obj, sizeof(__uint128_t) );
	m_ptr += sizeof(__uint128_t);
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack<__uint128_t>(unsigned char* m_ptr, __uint128_t &m_obj)
{
	memcpy( &m_obj, m_ptr, sizeof(__uint128_t) );
	m_ptr += sizeof(__uint128_t);
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for Date
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const Date &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		DATE_MEMBERS
	#undef VARIABLE

	return ret;	
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, Date &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
		DATE_MEMBERS
        #undef VARIABLE
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const Date &m_obj)
{

	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
		DATE_MEMBERS
	#undef VARIABLE

	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for SrirachaOptions
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const SrirachaOptions &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		SRIRACHA_OPTION_MEMBERS
	#undef VARIABLE

	return ret;	
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, SrirachaOptions &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
		SRIRACHA_OPTION_MEMBERS
        #undef VARIABLE
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const SrirachaOptions &m_obj)
{

	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
		SRIRACHA_OPTION_MEMBERS
	#undef VARIABLE

	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for SearchMatch
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const SearchMatch &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		SEARCH_MATCH_MEMBERS
	#undef VARIABLE

	return ret;	
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, SearchMatch &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
		SEARCH_MATCH_MEMBERS
        #undef VARIABLE
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const SearchMatch &m_obj)
{

	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
		SEARCH_MATCH_MEMBERS
	#undef VARIABLE

	return m_ptr;
}