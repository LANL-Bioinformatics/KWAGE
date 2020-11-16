#include "parse_tar.h"
#include <string.h>
#include <iostream>
#include <deque>

using namespace std;

size_t compute_file_size(unsigned char m_file_size[]);
size_t get_base256(unsigned char m_file_size[]);
size_t get_base8(unsigned char m_file_size[]);

TarIterator::TarIterator(const string &m_filename)
{
	// Use zlib to read both compressed and uncompressed tar files.
	fin = gzopen(m_filename.c_str(), "r");

	if(fin == NULL){

		cerr << "Error opening: " << m_filename << endl;
		throw __FILE__ ":TarIterator::TarIterator: Unable to open tar file";
	}

	bytes_to_next_header = 0;
	bytes_to_pad = 0;
	
	next();
}

void TarIterator::next()
{
	buffer.clear();

	if(bytes_to_next_header == 0){

		// Skip any padding bytes
		if(bytes_to_pad > 0){

			char *temp = new char [bytes_to_pad];

			if(temp == NULL){
				throw __FILE__ ":TarIterator: Unable to allocate temp buffer for skipping pad";
			}

			const int zret = gzread(fin, temp, bytes_to_pad);

			if( zret != int(bytes_to_pad) ){
				throw __FILE__ ":TarIterator::next: Error reading padding bases";
			}

			delete [] temp;
			temp = NULL;
		}

		// Read the tar file header
		const int zret = gzread( fin, &header, sizeof(TarHeader) );

		if(zret == 0){

			gzclose(fin);
			fin = NULL;

			return;
		}

		if( zret != int( sizeof(TarHeader) ) ){
			throw __FILE__ ":TarIterator::next: Error reading tar header";
		}

		// Is the header completely zeroed out?
		unsigned char *ptr = (unsigned char*)&header;

		bool all_zero = true;

		for(size_t i = 0;i < sizeof(TarHeader);++i){

			if(*ptr != 0){

				all_zero = false;
				break;
			}
		}

		// A completely empty header means it is time to stop reading
		if(all_zero){

			gzclose(fin);
			fin = NULL;

			return;
		}

		// Validate the header
		if( strncmp(header.magic, "ustar", 5) != 0 ){
			throw __FILE__ ":TarIterator::next: Error reading tar header magic value";
		}

		// It appears that some tar programs set the version to values other than "00"
		//if( strncmp(header.version, "00", 2) != 0 ){
		//	throw __FILE__ ":TarIterator::next: Error reading tar header version";
		//}

		// Compute the file size
		bytes_to_next_header = compute_file_size(header.file_size);
		bytes_to_pad = (512 - bytes_to_next_header%512)%512;

		// Convert the C-string into a C++ string
		file_name = string(header.file_name);
		
		// Attempt to load the first line of the first file. This has
		// the effect of skipping zero length entries (like directories).
		next();
	}
    else{
    
		deque<char> local;
		
		// A previous version of TarCat used gzgets(). However, this function can be derailed when
		// it encounters non-ASCII characters in a file (perhaps due to a corrupted file).

		// Limit the number of bytes read to the smaller of bytes_to_next_header + 1 (for the '\0')
		// and buffer_len. This required to handle files that are missing the '\n' on the last line.
		// Without a '\n', the gzgets() will read into the next record.
		//if(gzgets( fin, temp, min(bytes_to_next_header + 1, buffer_len) ) == NULL){
		//	throw __FILE__ ":TarIterator::next: Unexpected end of file";
		//}

		while(bytes_to_next_header > 0){
			
			const char c = gzgetc(fin);

			if(c == -1){
				throw __FILE__ ":TarIterator::next: Error reading byte!";
			}

			--bytes_to_next_header;

			if(c == '\n'){
				break;
			}

			// Do not save DOS carriage return symbols
			if(c != '\r'){
				local.push_back(c);
			}
		}

		buffer.assign( local.begin(), local.end() );
	}
}

size_t compute_file_size(unsigned char m_file_size[])
{
	// The high order bit of the left-most byte
	// to determine the encoding scheme
	if( (m_file_size[0] >> 7) & 1){
		return get_base256(m_file_size);
	}

	return get_base8(m_file_size);
}

size_t get_base256(unsigned char m_file_size[])
{
	size_t ret = 0;

	size_t place = 1;

	// Don't include m_file_size[0], as this byte
	// is not used (since it contained the bit that
	// indicated 256 byte encoding)
	for(int i = 11;i > 0;--i){

		ret += place*m_file_size[i];
		place *= 256;
	}

	return ret;
}

size_t get_base8(unsigned char m_file_size[])
{
	size_t ret = 0;

	if( (m_file_size[11] != ' ') && (m_file_size[11] != '\0') ){
		throw __FILE__ ":get_base8: Missing end of string symbol";
	}

	size_t place = 1;

	for(int i = 10;i >= 0;--i){

		switch(m_file_size[i]){
			case '0':
			break;
		case '1':
			ret += place*1;
			break;
		case '2':
			ret += place*2;
			break;
		case '3':
			ret += place*3;
			break;
		case '4':
			ret += place*4;
			break;
		case '5':
			ret += place*5;
			break;
		case '6':
			ret += place*6;
			break;
		case '7':
			ret += place*7;
			break;
		default:
			throw __FILE__ ":get_base8: Invalid symbol";
		};

		place *= 8;
	}

	return ret;
}
