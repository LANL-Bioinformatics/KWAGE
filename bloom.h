#ifndef __BLOOM_FILTER
#define __BLOOM_FILTER

#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <unordered_map>

//#include <iostream> // DEBUG

#include "hash.h"
#include "mpi_util.h"

#define		MAP				std::unordered_map
#define		MULTIMAP			std::unordered_multimap

#define		MIN_LOG_2_BLOOM_FILTER_LEN	20
#define		MAX_LOG_2_BLOOM_FILTER_LEN	31

#define		MIN_NUM_HASH			1
#define		MAX_NUM_HASH			5

// We need to protect any commas that appear in template variables
// if we are going to combine them with X Macros
#define SINGLE_ARG(...) __VA_ARGS__

class BitVector
{
	public:
		typedef unsigned char BLOCK;
		
		enum {
			BITS_PER_BLOCK = sizeof(BLOCK)*8,
			BLOCK_MASK = 0xFF
		};
		
	private:
				
		BLOCK* buffer;
		uint32_t num_bits;
	
	public:
		BitVector() : 
			buffer(NULL), num_bits(0)
		{
		};
		
		BitVector(const uint32_t &m_num_bits) : 
			num_bits(m_num_bits)
		{
			buffer = NULL;
			
			if(num_bits > 0){
				
				buffer = new BLOCK[ num_block() ];
				
				if(buffer == NULL){
					throw __FILE__ ":BitVector: Unable to allocate buffer";
				}
			}
		};
		
		~BitVector()
		{
			clear();
		}
		
		BitVector(const BitVector &m_rhs): buffer(NULL), num_bits(0)
		{
			*this = m_rhs;
		};
		
		inline void clear()
		{
			if(buffer != NULL){
			
				delete [] buffer;
				buffer = NULL;
			}
		};
		
		inline uint32_t get_num_bits() const
		{
			return num_bits;
		};
		
		inline uint32_t num_block() const
		{
			return num_bits/BITS_PER_BLOCK + ( (num_bits%BITS_PER_BLOCK) > 0 ? 1 : 0);
		};
		
		inline BLOCK* ptr()
		{
			// For combining distributed BitVectors with MPI Reduce
			return buffer;
		};
		
		// No bounds checking for fast access!
		inline bool get_bit(const size_t &m_index) const
		{
			// DEBUG
			//if(m_index >= num_bits){
			//	std::cerr << "m_index = " << m_index << std::endl;
			//	std::cerr << "num_bits = " << num_bits << std::endl;
			//	throw __FILE__ ":get_bit: Overflow!";
			//}
			
			// Within a block, bits are packed from high to low:
			// |76543210|76543210|76543210|76543210|76543210|76543210|...
			// |Block 0 |Block 1 |Block 2 |Block 3 |Block 4 |Block 5 |..
			return (buffer[m_index/BITS_PER_BLOCK] >> m_index%BITS_PER_BLOCK) & BLOCK(1);
		};
		
		// No bounds checking for fast access!
		inline void set_bit(const size_t &m_index)
		{
			// DEBUG
			//if(m_index/BITS_PER_BLOCK >= num_block()){
			//	
			//	std::cerr << "num_block() = " << num_block() << std::endl;
			//	std::cerr << "m_index = " << m_index << std::endl;
			//	std::cerr << "m_index/BITS_PER_BLOCK = " << m_index/BITS_PER_BLOCK << std::endl;
			//	std::cerr << "num_bits = " << num_bits << std::endl;
			//	throw __FILE__ ":set_bit: Overflow";
			//}

			// Within a block, bits are packed from high to low:
			// |76543210|76543210|76543210|76543210|76543210|76543210|...
			// |Block 0 |Block 1 |Block 2 |Block 3 |Block 4 |Block 5 |..
			buffer[m_index/BITS_PER_BLOCK] |= BLOCK(1) << m_index%BITS_PER_BLOCK;			
		};
		
		inline void unset_bit(const size_t &m_index)
		{
			// Within a block, bits are packed from high to low:
			// |76543210|76543210|76543210|76543210|76543210|76543210|...
			// |Block 0 |Block 1 |Block 2 |Block 3 |Block 4 |Block 5 |..
			buffer[m_index/BITS_PER_BLOCK] &= ~(BLOCK(1) << m_index%BITS_PER_BLOCK);
		};
		
		// Set *all* of the bits to zero
		inline void unset_all_bits()
		{
			if(buffer != NULL){
				memset( buffer, BLOCK(0), num_block()*sizeof(BLOCK) );
			}
		};
		
		// Set *all* of the bits to zero
		inline void set_all_bits()
		{
			if(buffer != NULL){
				memset( buffer, ~BLOCK(0), num_block()*sizeof(BLOCK) );
			}
		};

		// Enable deep copy construction
		inline BitVector& operator=(const BitVector &m_rhs)
		{
			if(buffer != NULL){
				delete [] buffer;
			}
			
			num_bits = m_rhs.num_bits;
			
			buffer = new BLOCK[ num_block() ];
				
			if(buffer == NULL){
				throw __FILE__ ":BitVector=(): Unable to allocate buffer";
			}
			
			memcpy( buffer, m_rhs.buffer, num_block() );
			
			return *this;
		};
		
		// Combine two-bit vectors with bitwise OR
		inline BitVector& operator|=(const BitVector &m_rhs)
		{
			if(num_bits != m_rhs.num_bits){
				
				// DEBUG
				//std::cerr << "buffer = " << size_t(buffer) << std::endl;
				//std::cerr << "num_bits = " << num_bits << std::endl;
				//std::cerr << "m_rhs.num_bits = " << m_rhs.num_bits << std::endl;

				if(buffer != NULL){
					delete [] buffer;
				}

				num_bits = m_rhs.num_bits;

				buffer = new BLOCK[ num_block() ];

				if(buffer == NULL){
					throw __FILE__ ":BitVector|=(): Unable to allocate buffer";
				}

				memcpy( buffer, m_rhs.buffer, num_block() );
				
				return *this;
			}
			
			// num_bits == m_rhs.num_bits
			for(uint32_t i = 0;i < num_block();++i){
				buffer[i] |= m_rhs.buffer[i];
			}
			
			return *this;
		};
		
		// Combine two-bit vectors with bitwise AND
		inline BitVector& operator&=(const BitVector &m_rhs)
		{
			if(num_bits != m_rhs.num_bits){
				
				// DEBUG
				//std::cerr << "num_bits = " << num_bits << std::endl;
				//std::cerr << " m_rhs.num_bits = " <<  m_rhs.num_bits << std::endl;
				
				throw __FILE__ ":BitVector::operator&=: Unequal lengths";
			}
			
			// num_bits == m_rhs.num_bits
			for(uint32_t i = 0;i < num_block();++i){
				buffer[i] &= m_rhs.buffer[i];
			}
			
			return *this;
		};
		
		inline void resize(const size_t &m_num_bits)
		{
			clear();
			
			num_bits = m_num_bits;
			
			buffer = new BLOCK[ num_block() ];
				
			if(buffer == NULL){
				throw __FILE__ ":BitVector::resize: Unable to allocate buffer";
			}
		};
		
		inline void increment_count(std::vector<unsigned int> &m_count) const
		{
			if(m_count.size() != num_bits){
				m_count.resize(num_bits, 0);
			}
			
			const uint32_t last_block = num_block() - 1;
			
			// Within a block, bits are packed from high to low:
			// |76543210|76543210|76543210|76543210|76543210|76543210|...
			// |Block 0 |Block 1 |Block 2 |Block 3 |Block 4 |Block 5 |..
			for(uint32_t i = 0;i < last_block;++i){
				
				// We only need to increment the count of a block has
				// one or more non-zero bits
				if(buffer[i]){
					
					std::vector<unsigned int>::iterator iter = 
						m_count.begin() + i*BITS_PER_BLOCK;
					
					for(uint32_t j = 0;j < BITS_PER_BLOCK;++j, ++iter){
						*iter += (buffer[i] >> j) & 1;
					}
				}
			}
			
			// Handle the last block carefully since not every bit is guaranteed to be valid
			if(buffer[last_block]){
					
				std::vector<unsigned int>::iterator iter = 
					m_count.begin() + last_block*BITS_PER_BLOCK;
				
				// Stop checking bits when we have tested the last valid element in m_count
				for(uint32_t j = 0;(j < BITS_PER_BLOCK) && ( iter != m_count.end() );
					++j, ++iter){
					
					*iter += (buffer[last_block] >> j) & 1;
				}
			}
		};
		
		// Is any bit set in the bit vector filter?
		inline bool max_bit()
		{
			const size_t last_block = num_block() - 1;
			
			for(size_t i = 0;i < last_block;++i){
				
				// We only need to increment the count of a block has
				// one or more non-zero bits
				if(buffer[i]){
					return true;
				}
			}
			
			// Handle the last block
			const unsigned char remainder = num_bits%BITS_PER_BLOCK;
			
			if( (remainder == 0) && buffer[last_block] ){
				return true;
			}
			
			for(unsigned char i = 0;i < remainder;++i){
				if( (buffer[last_block] >> i) & 1 ){
					return true;
				}
			}
			
			return false;
		};
		
		inline void assign_bits(const uint32_t &m_num_bits, const BitVector &m_src)
		{
			clear();
			
			num_bits = m_num_bits;
			
			buffer = new BLOCK[ num_block() ];
				
			if(buffer == NULL){
				throw __FILE__ ":BitVector::assign_bits: Unable to allocate buffer";
			}
			
			if(num_bits > m_src.num_bits){
				throw __FILE__ ":BitVector::assign_bits: |dst| > |src|";
			}
			
			const uint32_t last_block =  num_block() - 1;
			
			memcpy( buffer, m_src.buffer, last_block);
			
			// Handle the last block
			const unsigned char remainder = num_bits%BITS_PER_BLOCK;
			
			if(remainder == 0){
				buffer[last_block] = m_src.buffer[last_block];
			}
			else{
				buffer[last_block] = 0;
				
				for(unsigned char i = 0;i < remainder;++i){
					if( (m_src.buffer[last_block] >> i) & 1 ){
						buffer[last_block] |= 1 << i;
					}
				}
			}
		};
		
		// Count the number of bits that have been set to 1
		inline size_t count() const
		{
			
			size_t ret = 0;
			
			// For efficient bit counting, use the __builtin_popcountl intrinsic function
			const uint32_t len = ( (num_block() - 1)*sizeof(BLOCK) )/sizeof(unsigned long int);
			
			unsigned long int *ptr = (unsigned long int*)buffer;
			
			for(uint32_t i = 0;i < len;++i, ++ptr){
				ret += __builtin_popcountl(*ptr);
			}
			
			const uint32_t remainder_bits_begin = len*sizeof(unsigned long int)*8;
			
			for(uint32_t i = remainder_bits_begin;i < num_bits;++i){
				ret += (get_bit(i) ? 1 : 0);
			}
			
			return ret;
		};
		
		inline bool empty() const
		{
			return (num_bits == 0);
		};
		
		unsigned int crc32() const;
		
		// Write to a file
		inline void write(std::ofstream &m_fout) const
		{
			m_fout.write( (char*)buffer, num_block()*sizeof(BLOCK) );
		};
		
		// Write to a memory buffer
		inline void write(unsigned char *m_ptr) const
		{
			memcpy(m_ptr, buffer, num_block()*sizeof(BLOCK) );
		};
		
		// Read from a file
		inline void read(std::ifstream &m_fin)
		{
			m_fin.read( (char*)buffer, num_block()*sizeof(BLOCK) );
		};
		
		// Read from a memory buffer
		inline void read(unsigned char *m_ptr)
		{
			memcpy(buffer, m_ptr, num_block()*sizeof(BLOCK) );
		};

		template<class T> friend size_t mpi_size(const T &m_obj);
		template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr,
			const T &m_obj);
		template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, 
			T &m_obj);

		template<class T> friend void binary_write(std::ostream &m_out,
			const T &m_obj);
		template<class T> friend void binary_read(std::istream &m_in, 
			T &m_obj);
};

template<> size_t mpi_size(const BitVector &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const BitVector &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, BitVector &m_obj);

template<> void binary_write(std::ostream &m_out, const BitVector &m_obj);
template<> void binary_read(std::istream &m_in, BitVector &m_obj);

// Store all of the SRA metadata associated with a Bloom filter. Use the mpi pack and unpack
// friend functions for disk-based serialization
struct FilterInfo
{
	// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
	// ensure that structure variables are correctly serialized
	#define FILTER_INFO_MEMBERS \
		VARIABLE(std::string, run_accession) \
		VARIABLE(std::string, experiment_accession) \
		VARIABLE(std::string, experiment_title) \
		VARIABLE(std::string, experiment_design_description) \
		VARIABLE(std::string, experiment_library_name) \
		VARIABLE(std::string, experiment_library_strategy) \
		VARIABLE(std::string, experiment_library_source) \
		VARIABLE(std::string, experiment_library_selection) \
		VARIABLE(std::string, experiment_instrument_model) \
		VARIABLE(std::string, sample_accession) \
		VARIABLE(std::string, sample_taxa) \
		VARIABLE(SINGLE_ARG(MULTIMAP<std::string, std::string>), sample_attributes) \
		VARIABLE(std::string, study_accession) \
		VARIABLE(std::string, study_title) \
		VARIABLE(std::string, study_abstract)
	
	#define VARIABLE(A, B) A B;
		FILTER_INFO_MEMBERS
	#undef VARIABLE

	std::string csv_string() const;
	std::string json_string(const std::string &m_prefix) const;

	template<class T> friend size_t mpi_size(const T &m_obj);
	template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr,
		const T &m_obj);
	template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, 
		T &m_obj);

	template<class T> friend void binary_write(std::ostream &m_out,
		const T &m_obj);
	template<class T> friend void binary_read(std::istream &m_in, 
		T &m_obj);
};

template<> size_t mpi_size(const FilterInfo &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const FilterInfo &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, FilterInfo &m_obj);

template<> void binary_write(std::ostream &m_out, const FilterInfo &m_obj);
template<> void binary_read(std::istream &m_in, FilterInfo &m_obj);

struct BloomParam
{
	// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
	// ensure that structure variables are correctly serialized
	#define BLOOM_PARAM_MEMBERS \
		VARIABLE(uint32_t, kmer_len) \
		VARIABLE(uint32_t, log_2_filter_len) \
		VARIABLE(uint32_t, num_hash) \
		VARIABLE(HashFunction, hash_func) \
	
	#define VARIABLE(A, B) A B;
		BLOOM_PARAM_MEMBERS
	#undef VARIABLE

	BloomParam()
	{
		#define VARIABLE(A, B) B = 0;
			BLOOM_PARAM_MEMBERS
		#undef VARIABLE
	};

	// Needed for using BloomParam as a hash function key
	inline bool operator==(const BloomParam &m_rhs) const
	{
		return (kmer_len == m_rhs.kmer_len) &&
			(log_2_filter_len == m_rhs.log_2_filter_len) && 
			(num_hash == m_rhs.num_hash) && 
			(hash_func == m_rhs. hash_func);
	};

	inline bool operator<(const BloomParam &m_rhs) const
	{
		if(kmer_len == m_rhs.kmer_len){

			if(log_2_filter_len == m_rhs.log_2_filter_len){
				return (num_hash < m_rhs.num_hash);
			}

			return (log_2_filter_len < m_rhs.log_2_filter_len);
		}

		return kmer_len < m_rhs.kmer_len;
	};

	// Return a size_t to allow for very large Bloom filters
	inline size_t filter_len() const
	{
		return (size_t(1) << log_2_filter_len);
	};
	
	inline bool empty() const
	{
		return (kmer_len == 0) &&
			(log_2_filter_len == 0) &&
			(num_hash == 0);
	};
	
	template<class T> friend size_t mpi_size(const T &m_obj);
	template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr,
		const T &m_obj);
	template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, 
			T &m_obj);

	template<class T> friend void binary_write(std::ostream &m_out,
		const T &m_obj);
	template<class T> friend void binary_read(std::istream &m_in, 
			T &m_obj);
};

template<> size_t mpi_size(const BloomParam &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const BloomParam &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, BloomParam &m_obj);

template<> void binary_write(std::ostream &m_out, const BloomParam &m_obj);
template<> void binary_read(std::istream &m_in, BloomParam &m_obj);

// Allow the use of the log_2_filter_len and num_hash values from the BloomParam
// structure to serve as a hash key. For now, we do not include a contribution
// from the kmer length.
namespace std
{
	template<>
	struct hash <BloomParam>
	{
	    size_t operator()(const BloomParam& m_param) const
	    {
        	return (size_t(m_param.num_hash) << 32) | size_t(m_param.log_2_filter_len);
	    }
	};
}

class BloomFilter : public BitVector
{
	private:
		
		BloomParam param;
		unsigned int bloom_crc32;
		FilterInfo info;
		
	public:
		BloomFilter() : bloom_crc32(0)
		{
		};
		
		BloomFilter(const BloomParam &m_param) : 
			BitVector( m_param.filter_len() ),
			param(m_param), bloom_crc32(0)
		{
		};
		
		~BloomFilter()
		{
			clear();
		}
		
		BloomFilter(const BloomFilter &m_rhs)
		{
			*this = m_rhs;
		};
		
		inline void set_info(const FilterInfo &m_info)
		{
			info = m_info;
		};

		inline const FilterInfo& get_info() const
		{
			return info;
		};
		
		inline void set_param(const BloomParam &m_param)
		{
			if(param.log_2_filter_len != m_param.log_2_filter_len){

				// Resize the BitVector to be consistent with the new
				// Bloom filter size and unset all bits (so we can accumulate
				// Bloom filters with binary OR)
				BitVector::resize( m_param.filter_len() );
				BitVector::unset_all_bits();
			}

			param = m_param;
		};

		inline const BloomParam& get_param() const
		{
			return param;
		};
		
		// Enable deep copy construction
		inline BloomFilter& operator=(const BloomFilter &m_rhs)
		{
			param = m_rhs.param;
			info = m_rhs.info;
			
			BitVector::operator=(m_rhs);
			
			return *this;
		};
		
		unsigned int update_crc32();
		bool test_crc32() const;
		
		inline unsigned int get_crc32() const
		{
			return bloom_crc32;
		};
		
		template<class T> friend size_t mpi_size(const T &m_obj);
		template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr,
			const T &m_obj);
		template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, 
			T &m_obj);

		template<class T> friend void binary_write(std::ostream &m_out,
			const T &m_obj);
		template<class T> friend void binary_read(std::istream &m_in, 
			T &m_obj);
};

template<> size_t mpi_size(const BloomFilter &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const BloomFilter &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, BloomFilter &m_obj);

template<> void binary_write(std::ostream &m_out, const BloomFilter &m_obj);
template<> void binary_read(std::istream &m_in, BloomFilter &m_obj);

class SubFilter : public BitVector
{
	public:
		enum {
			INVALID_FILE = 0xFFFFFFFF,
			INVALID_LOC = 0xFFFFFFFF
		};
		
	private:
	
		unsigned int file_index;
		unsigned int src_loc; // Either file or memory location
		unsigned int dst_loc; // Order for output

	public:
		SubFilter() : file_index(INVALID_FILE), src_loc(INVALID_LOC), dst_loc(INVALID_LOC)
		{
		};
		
		~SubFilter()
		{
		}
		
		SubFilter(const SubFilter &m_rhs)
		{
			*this = m_rhs;
		};
				
		// Enable deep copy construction
		inline SubFilter& operator=(const SubFilter &m_rhs)
		{
			file_index = m_rhs.file_index;
			src_loc = m_rhs.src_loc;
			dst_loc = m_rhs.dst_loc;
			
			BitVector::operator=(m_rhs);
			
			return *this;
		};

		// Enable the sorting of SubFilters based on their
		// destination index
		inline bool operator<(const SubFilter &m_rhs) const
		{
			return dst_loc < m_rhs.dst_loc;
		};

		inline void set_src_loc(const unsigned int &m_loc)
		{
			src_loc = m_loc;
		};
		
		inline unsigned int get_src_loc() const
		{
			return src_loc;
		};
		
		inline void set_dst_loc(const unsigned int &m_loc)
		{
			dst_loc = m_loc;
		};
		
		inline unsigned int get_dst_loc() const
		{
			return dst_loc;
		};
		
		inline void set_file(const unsigned int &m_index)
		{
			file_index = m_index;
		};
		
		inline unsigned int get_file() const
		{
			return file_index;
		};

		template<class T> friend size_t mpi_size(const T &m_obj);
		template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr,
			const T &m_obj);
		template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, 
			T &m_obj);
};

template<> size_t mpi_size(const SubFilter &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const SubFilter &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, SubFilter &m_obj);

BloomParam optimal_bloom_param(const uint32_t &m_kmer_len, const size_t &m_num_kmer, 
	const float &m_p, const HashFunction &m_func);

#endif // __BLOOM_FILTER
