#ifndef __CALDERA
#define __CALDERA

#define	CALDERA_VERSION		"0.4c"
#define	INVENTORY_VERSION	"0.7"
#define	MAESTRO_VERSION		"0.9"

#include <cstdint> // uint32_t
#include <fstream>
#include "hash.h"

// Needed to supress compiler wanings on g++
#include "binary_io.h"

// The currently allowed compression levels
enum {
	NO_COMPRESSION,
	RLE_COMPRESSION,
	RLE_HUFFMAN_COMPRESSION
};

#define	CALDERA_MAGIC_NUMBER			0x20191025
#define	CURRENT_DBFILE_VERSION		2

// MPI Messages
#define	INFO_BUFFER			1001
#define	INFO_LEN			1002

// The file header for the CALDERA database files
struct DBFileHeader
{
	// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
	// ensure that structure variable are correctly serialized
	#define DBFILEHEADER_MEMBERS \
		VARIABLE(uint32_t, magic) \
		VARIABLE(uint32_t, version) \
		VARIABLE(uint32_t, crc32) \
		VARIABLE(uint32_t, kmer_len) \
		VARIABLE(uint32_t, num_hash) \
		VARIABLE(uint32_t, log_2_filter_len) \
		VARIABLE(uint32_t, num_filter) \
		VARIABLE(HashFunction, hash_func) \
		VARIABLE(uint32_t, compression) \
		VARIABLE(uint64_t, info_start)

	#define VARIABLE(A, B) A B;
		DBFILEHEADER_MEMBERS
	#undef VARIABLE
	
	DBFileHeader()
	{
		// Initialize all variables to zero (which will work until
		// we add a member variable that can not be assigned to zero).
		#define VARIABLE(A, B) B = 0;
			DBFILEHEADER_MEMBERS
		#undef VARIABLE

		magic = CALDERA_MAGIC_NUMBER;
		version = CURRENT_DBFILE_VERSION;
	};
	
	// Return a size_t to allow for very large Bloom filters
	inline size_t filter_len() const
	{
		return (size_t(1) << log_2_filter_len);
	};
	
	template<class T> friend void binary_write(std::ostream &m_out,
			const T &m_obj);
	template<class T> friend void binary_read(std::istream &m_in, 
		T &m_obj);
};

template<> void binary_write(std::ostream &m_out, const DBFileHeader &m_obj);
template<> void binary_read(std::istream &m_in, DBFileHeader &m_obj);

#endif // __CALDERA
