// Classes for managing slice compression (deflation) and slice decompress (inflation)
#ifndef __SLICE_Z
#define __SLICE_Z

#include <fstream>
#include <zlib.h>

#define		MAX_COMPRESSED_BYTES	256 // (i.e. NUM_GROUP_CHUNK/( 8 bits per byte) )
#define		ZLIB_WINDOW_BITS	-9 // No zlib or gzip header

template<unsigned int MAX_BUFFER_SIZE>
class InflateSlice
{
    private:
	z_stream engine;
	unsigned char *z_src_buffer;
	unsigned char *z_dst_buffer;
	unsigned int z_src_len;
	unsigned int z_dst_len;

    public:

        InflateSlice()
        {
		engine.zalloc = Z_NULL;
		engine.zfree = Z_NULL;
		engine.opaque = Z_NULL;
		engine.avail_in = 0;
		engine.next_in = NULL;
		engine.avail_out = MAX_BUFFER_SIZE;
		engine.next_out = NULL;
		
		z_src_buffer = new unsigned char[MAX_BUFFER_SIZE];

		if(z_src_buffer == NULL){
			throw __FILE__ ":InflateSlice: Unable to allocate z_src_buffer";
		};

		z_dst_buffer = new unsigned char[MAX_BUFFER_SIZE];

		if(z_dst_buffer == NULL){
			throw __FILE__ ":InflateSlice: Unable to allocate z_src_buffer";
		};

		if( Z_OK != inflateInit2(&engine, ZLIB_WINDOW_BITS) ){
			throw __FILE__ ":InflateSlice: Error calling inflateInit2";
		}
        };

	InflateSlice(const InflateSlice &m_rhs)
	{
		z_src_buffer = new unsigned char[MAX_BUFFER_SIZE];

		if(z_src_buffer == NULL){
			throw __FILE__ ":InflateSlice: Unable to allocate z_src_buffer";
		};

		z_dst_buffer = new unsigned char[MAX_BUFFER_SIZE];

		if(z_dst_buffer == NULL){
			throw __FILE__ ":InflateSlice: Unable to allocate z_src_buffer";
		};
		
		*this = m_rhs;
	};
	
        ~InflateSlice()
        {
		if(z_src_buffer != NULL){

			delete [] z_src_buffer;
			z_src_buffer = NULL;
		}

		if(z_dst_buffer != NULL){

			delete [] z_dst_buffer;
			z_dst_buffer = NULL;
		}

		if( Z_OK != inflateEnd(&engine) ){
			throw __FILE__ ":~InflateSlice: Error in inflateEnd";
		}
        };

        void inflate(std::ifstream &m_fin, unsigned char m_len)
        {
		m_fin.read( (char*)z_src_buffer, m_len);

		if(!m_fin){
			throw __FILE__":InflateSlice::inflate(fstream): Error reading compressed slice from file";
		}

		engine.next_in = z_src_buffer;
		engine.avail_in = m_len;

		engine.next_out = z_dst_buffer;
		engine.avail_out = MAX_COMPRESSED_BYTES;

		if( Z_STREAM_END != ::inflate(&engine, Z_FINISH) ){
			throw __FILE__ ":InflateSlice::inflate(fstream): Error in inflate (1)";
		}

		if( Z_OK != inflateReset(&engine) ){
			throw __FILE__ ":InflateSlice::inflate(fstream): Error in inflateReset (1)";
		}
        };

        void inflate(unsigned char* m_buffer, unsigned char m_len)
        {
		engine.next_in = m_buffer;
		engine.avail_in = m_len;

		engine.next_out = z_dst_buffer;
		engine.avail_out = MAX_COMPRESSED_BYTES;

		if( Z_STREAM_END != ::inflate(&engine, Z_FINISH) ){
			throw __FILE__ ":InflateSlice::inflate(unsigned char*): Error in inflate (1)";
		}

		if( Z_OK != inflateReset(&engine) ){
			throw __FILE__ ":InflateSlice::inflate(unsigned char*): Error in inflateReset (1)";
		}
        };

        unsigned char* ptr() const
        {
            return z_dst_buffer;
        };

        inline unsigned int size() const
        {
            return MAX_BUFFER_SIZE - engine.avail_out;
        };
	
	// Do not allow Inflate slice to be copied (copying will require 
	// some work to manage the decompression engine state)!
	InflateSlice& operator=(const InflateSlice &m_rhs)
	{
		memcpy( z_src_buffer, m_rhs.z_src_buffer, MAX_BUFFER_SIZE*sizeof(unsigned char) );
		memcpy( z_dst_buffer, m_rhs.z_dst_buffer, MAX_BUFFER_SIZE*sizeof(unsigned char) );
		
		z_src_len = m_rhs.z_src_len;
		z_dst_len = m_rhs.z_dst_len;
		
		inflateCopy(&engine, (z_stream*)&m_rhs.engine);
		
		return *this;
	};
};

template<unsigned int MAX_BUFFER_SIZE>
class CompressSlice
{
    private:
        z_stream engine;
        unsigned char *z_buffer;
        unsigned int z_buffer_len;

    public:

        CompressSlice(const int &m_level = 6, const int &m_memLevel = 9)
        {
            engine.zalloc = Z_NULL;
            engine.zfree = Z_NULL;
            engine.opaque = Z_NULL;
            engine.avail_in = 0;
            engine.next_in = NULL;
            engine.avail_out = MAX_BUFFER_SIZE;
	    engine.next_out = NULL;

            z_buffer = new unsigned char[MAX_BUFFER_SIZE];

            if(z_buffer == NULL){
                throw __FILE__ ":CompressSlice: Unable to allocate z_buffer";
            };

            if( Z_OK != deflateInit2(&engine, 
                m_level,
                Z_DEFLATED,
                ZLIB_WINDOW_BITS, 
                m_memLevel, 
                Z_RLE) ){ // Only Z_RLE is supported for now

                throw __FILE__ ":CompressSlice: Error calling deflateInit2";
            }
        };

	CompressSlice(const CompressSlice &m_rhs)
	{
		z_buffer = new unsigned char[MAX_BUFFER_SIZE];

		if(z_buffer == NULL){
			throw __FILE__ ":CompressSlice: Unable to allocate z_buffer";
		};
		
		*this = m_rhs;
	};
	
        ~CompressSlice()
        {
            if(z_buffer != NULL){

                delete [] z_buffer;
                z_buffer = NULL;
            }

            // Don't check the return value of deflateEnd, since
            // we may have discarded some output when choosing not
            // to use compression.
            deflateEnd(&engine);
        };

        unsigned char* ptr() const
        {
            return z_buffer;
        };

        inline unsigned int size() const
        {
            return MAX_BUFFER_SIZE - engine.avail_out;
        };

        inline bool compress(unsigned char* m_ptr, unsigned char m_len)
        {
            engine.avail_in = m_len;
            engine.next_in = m_ptr;

            engine.avail_out = MAX_BUFFER_SIZE;
            engine.next_out = z_buffer;

            const int zlib_ret = deflate(&engine, Z_FINISH);

            if( Z_OK != deflateReset(&engine) ){
                throw __FILE__ ":CompressSlice::compress: Error in deflateReset";
            }

            // Were we able to compresss?
            return ( (Z_STREAM_END == zlib_ret) && (size() < m_len) );
        };
	
	CompressSlice& operator=(const CompressSlice &m_rhs)
	{
		memcpy( z_buffer, m_rhs.z_buffer, MAX_BUFFER_SIZE*sizeof(unsigned char) );
		
		z_buffer_len = m_rhs.z_buffer_len;
		
		deflateCopy( &engine, (z_stream*)&(m_rhs.engine) );
		
		return *this;
	};
};

#endif // __SLICE_Z
