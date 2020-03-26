#include "binary_io.h"
#include "options.h"
#include "bloom.h"
#include "bigsi++.h"
#include "date.h"
#include <string.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for std::string
/////////////////////////////////////////////////////////////////////////////////////////
template<>
void binary_write<string>(ostream &m_out, const string &m_str)
{
	binary_write( m_out, m_str.size() );
	
	m_out.write( m_str.c_str(), m_str.size() );
	
	if(!m_out){
		throw __FILE__ ":binary_write<string>: Unable to write string";
	}
}

template<>
void binary_read<string>(istream &m_in, string &m_str)
{
	size_t len;
	
	binary_read(m_in, len);
	
	char* buffer = new char [len];

	if(buffer == NULL){
		throw __FILE__ ":binary_read<string>: Unable to allocate buffer";
	}

	m_in.read(buffer, len);

	if(!m_in){
		throw __FILE__ ":binary_read<string>: Unable to read string";
	}
	
	m_str.assign(buffer, len);
	
	delete [] buffer;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for __uint128_t
/////////////////////////////////////////////////////////////////////////////////////////

template<>
void binary_write<__uint128_t>(ostream &m_out, const __uint128_t &m_obj)
{
	m_out.write( (char*)&m_obj, sizeof(__uint128_t) );
	
	if(!m_out){
		throw __FILE__ ":binary_write<__uint128_t>: Unable to write __uint128_t";
	}
}

template<>
void binary_read<__uint128_t>(istream &m_in, __uint128_t &m_obj)
{
	m_in.read( (char*)&m_obj, sizeof(__uint128_t) );
	
	if(!m_in){
		throw __FILE__ ":binary_read<__uint128_t>: Unable to read __uint128_t";
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for BitVector
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
void binary_write(ostream &m_out, const BitVector &m_obj)
{
	binary_write(m_out, m_obj.num_bits);
	
	m_out.write( (char*)m_obj.buffer, m_obj.num_block() );
	
	if(!m_out){
		throw __FILE__ ":binary_write<BitVector>: Unable to write BitVector";
	}
}

template<> 
void binary_read(istream &m_in, BitVector &m_obj)
{
	m_obj.clear();
	
	binary_read(m_in, m_obj.num_bits);

	m_obj.resize(m_obj.num_bits);
	
	m_in.read( (char*)m_obj.buffer, m_obj.num_block() );
	
	if(!m_in){
		throw __FILE__ ":binary_read<BitVector>: Unable to read BitVector";
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for BloomParam
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
void binary_write(ostream &m_out, const BloomParam &m_obj)
{

	#define VARIABLE(A, B) binary_write(m_out, m_obj.B);
		BLOOM_PARAM_MEMBERS
	#undef VARIABLE
	
	if(!m_out){
		throw __FILE__ ":binary_write<BloomParam>: Unable to write BloomParam";
	}
}

template<> 
void binary_read(istream &m_in, BloomParam &m_obj)
{
	#define VARIABLE(A, B) binary_read(m_in, m_obj.B);
		BLOOM_PARAM_MEMBERS
	#undef VARIABLE
	
	if(!m_in){
		throw __FILE__ ":binary_read<BloomParam>: Unable to read BloomParam";
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for FilterInfo
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
void binary_write(ostream &m_out, const FilterInfo &m_obj)
{
	#define VARIABLE(A, B) binary_write(m_out, m_obj.B);
		FILTER_INFO_MEMBERS
	#undef VARIABLE
	
	if(!m_out){
		throw __FILE__ ":binary_write<FilterInfo>: Unable to write FilterInfo";
	}
}

template<> 
void binary_read(istream &m_in, FilterInfo &m_obj)
{
	#define VARIABLE(A, B) binary_read(m_in, m_obj.B);
		FILTER_INFO_MEMBERS
    	#undef VARIABLE
	
	if(!m_in){
		throw __FILE__ ":binary_read<FilterInfo>: Unable to read FilterInfo";
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for BloomFilter
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
void binary_write(ostream &m_out, const BloomFilter &m_obj)
{
	binary_write(m_out, m_obj.param);
	binary_write(m_out, m_obj.bloom_crc32);
	binary_write(m_out, m_obj.info);
	binary_write<BitVector>(m_out, m_obj);
	
	if(!m_out){
		throw __FILE__ ":binary_write<BloomFilter>: Unable to write BloomFilter";
	}
}

template<> 
void binary_read(istream &m_in, BloomFilter &m_obj)
{
	binary_read(m_in, m_obj.param);
	binary_read(m_in, m_obj.bloom_crc32);
	binary_read(m_in, m_obj.info);
	binary_read<BitVector>(m_in, m_obj);
	
	if(!m_in){
		throw __FILE__ ":binary_read<BloomFilter>: Unable to read BloomFilter";
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for DBFileHeader
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
void binary_write(ostream &m_out, const DBFileHeader &m_obj)
{
	#define VARIABLE(A, B) binary_write(m_out, m_obj.B);
		DBFILEHEADER_MEMBERS
	#undef VARIABLE
	
	if(!m_out){
		throw __FILE__ ":binary_write<DBFileHeader>: Unable to write DBFileHeader";
	}
}

template<> 
void binary_read(istream &m_in, DBFileHeader &m_obj)
{
	#define VARIABLE(A, B) binary_read(m_in, m_obj.B);
		DBFILEHEADER_MEMBERS
    	#undef VARIABLE
	
	if(!m_in){
		throw __FILE__ ":binary_read<DBFileHeader>: Unable to read DBFileHeader";
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for Date
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
void binary_write(ostream &m_out, const Date &m_obj)
{

	#define VARIABLE(A, B) binary_write(m_out, m_obj.B);
		DATE_MEMBERS
	#undef VARIABLE
	
	if(!m_out){
		throw __FILE__ ":binary_write<Date>: Unable to write Date";
	}
}

template<> 
void binary_read(istream &m_in, Date &m_obj)
{
	#define VARIABLE(A, B) binary_read(m_in, m_obj.B);
		DATE_MEMBERS
	#undef VARIABLE
	
	if(!m_in){
		throw __FILE__ ":binary_read<Date>: Unable to read Date";
	}
}