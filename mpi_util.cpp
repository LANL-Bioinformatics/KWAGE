#include "mpi_util.h"
#include "options.h"
#include "bloom.h"
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
// Specialization for BuildOptions
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const BuildOptions &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		BUILD_OPTION_MEMBERS
	#undef VARIABLE

	return ret;	
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, BuildOptions &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
		BUILD_OPTION_MEMBERS
        #undef VARIABLE
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const BuildOptions &m_obj)
{

	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
		BUILD_OPTION_MEMBERS
	#undef VARIABLE

	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for DownloadOptions
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const DownloadOptions &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		DOWNLOAD_OPTION_MEMBERS
	#undef VARIABLE

	return ret;	
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, DownloadOptions &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
		DOWNLOAD_OPTION_MEMBERS
        #undef VARIABLE
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const DownloadOptions &m_obj)
{

	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
		DOWNLOAD_OPTION_MEMBERS
	#undef VARIABLE

	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for BloomerOptions
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const BloomerOptions &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		BLOOMER_OPTION_MEMBERS
	#undef VARIABLE

	return ret;	
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, BloomerOptions &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
		BLOOMER_OPTION_MEMBERS
        #undef VARIABLE
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const BloomerOptions &m_obj)
{

	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
		BLOOMER_OPTION_MEMBERS
	#undef VARIABLE

	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for BitVector
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const BitVector &m_obj)
{
	return mpi_size(m_obj.num_bits) + m_obj.num_block();
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, BitVector &m_obj)
{
	m_obj.clear();
	
	m_ptr = mpi_unpack(m_ptr, m_obj.num_bits);

	m_obj.resize(m_obj.num_bits);
	
	memcpy( m_obj.buffer, m_ptr, m_obj.num_block() );
	m_ptr += m_obj.num_block();
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const BitVector &m_obj)
{
	m_ptr = mpi_pack(m_ptr, m_obj.num_bits);
	
	memcpy( m_ptr, m_obj.buffer, m_obj.num_block() );
	m_ptr += m_obj.num_block();
	
	return m_ptr;
}


/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for BloomFilter
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const BloomFilter &m_obj)
{
	return  mpi_size(m_obj.param) + 
		mpi_size(m_obj.bloom_crc32) + 
		mpi_size(m_obj.info) + 
		mpi_size<BitVector>(m_obj);
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, BloomFilter &m_obj)
{
	m_ptr = mpi_unpack(m_ptr, m_obj.param);
	m_ptr = mpi_unpack(m_ptr, m_obj.bloom_crc32);
	m_ptr = mpi_unpack(m_ptr, m_obj.info);
	m_ptr = mpi_unpack<BitVector>(m_ptr, m_obj);
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const BloomFilter &m_obj)
{
	m_ptr = mpi_pack(m_ptr, m_obj.param);
	m_ptr = mpi_pack(m_ptr, m_obj.bloom_crc32);
	m_ptr = mpi_pack(m_ptr, m_obj.info);
	m_ptr = mpi_pack<BitVector>(m_ptr, m_obj);
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for SubFilter
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const SubFilter &m_obj)
{
	return mpi_size(m_obj.file_index) + mpi_size(m_obj.src_loc) + mpi_size(m_obj.dst_loc) + 
		mpi_size<BitVector>(m_obj);
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, SubFilter &m_obj)
{
	m_ptr = mpi_unpack(m_ptr, m_obj.file_index);
	m_ptr = mpi_unpack(m_ptr, m_obj.src_loc);
	m_ptr = mpi_unpack(m_ptr, m_obj.dst_loc);
	m_ptr = mpi_unpack<BitVector>(m_ptr, m_obj);
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const SubFilter &m_obj)
{
	m_ptr = mpi_pack(m_ptr, m_obj.file_index);
	m_ptr = mpi_pack(m_ptr, m_obj.src_loc);
	m_ptr = mpi_pack(m_ptr, m_obj.dst_loc);
	m_ptr = mpi_pack<BitVector>(m_ptr, m_obj);
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for BloomParam
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const BloomParam &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		BLOOM_PARAM_MEMBERS
	#undef VARIABLE

	return ret;	
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, BloomParam &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
		BLOOM_PARAM_MEMBERS
        #undef VARIABLE
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const BloomParam &m_obj)
{

	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
		BLOOM_PARAM_MEMBERS
	#undef VARIABLE

	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for FilterInfo
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const FilterInfo &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		FILTER_INFO_MEMBERS
	#undef VARIABLE

	return ret;	
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, FilterInfo &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
		FILTER_INFO_MEMBERS
        #undef VARIABLE
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const FilterInfo &m_obj)
{

	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
		FILTER_INFO_MEMBERS
	#undef VARIABLE

	return m_ptr;
}
