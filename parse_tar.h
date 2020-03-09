#ifndef __PARSE_TAR
#define __PARSE_TAR

#include <string>
#include <zlib.h>

class TarIterator
{
	private:
		gzFile fin;
		std::string buffer;
		std::string file_name;
		size_t bytes_to_next_header; 
		size_t bytes_to_pad;

		// USTar format from Wikipedia: https://en.wikipedia.org/wiki/Tar_(computing)
		struct TarHeader
		{
			char file_name[100];
			char file_mode[8];
			char owner_id[8];
			char group_id[8];
			unsigned char file_size[12];
			char modification_time[12];
			char checksum[8];
			char type_flag;
			char linked_filename[100];
			char magic[6]; // ustar\0
			char version[2];
			char owner_user_name[32];
			char owner_group_name[32];
			char device_major_number[8];
			char device_minor_number[8];
			char file_name_prefix[155];
			char unused[12];
		} header;

		void next();
	public:

		TarIterator()
		{
			fin = NULL;
			bytes_to_next_header = 0;
			bytes_to_pad = 0;
		};

		TarIterator(const std::string &m_filename);

		~TarIterator()
		{
			if(fin != NULL){

				gzclose(fin);
				fin = NULL;
			}
		};

		inline operator bool() const
		{
			return (fin != NULL);
		};

		// Read the next line
		inline void operator++()
		{
			next();
		};

		inline const std::string& operator*() const
		{
			return buffer;
		};

		inline const std::string& filename() const
		{
			// We now keep the file name an internal std::string
			// to avoid the cost of converting from a C-string.
			//return std::string(header.file_name);
			return file_name;
		};
};

#endif // __PARSE_TAR
