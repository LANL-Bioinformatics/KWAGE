#include <fstream>
#include <iostream>
#include <deque>
#include <string.h>

#include "parse_sequence.h"
#include "file_util.h"

using namespace std;

FileType get_file_type(const string &m_filename);

FileType get_file_type(const string &m_filename)
{
	if( find_file_extension(m_filename, ".fna") || find_file_extension(m_filename, ".fna.gz") ||
	    find_file_extension(m_filename, ".fa") || find_file_extension(m_filename, ".fa.gz") ||
	    find_file_extension(m_filename, ".fasta") || find_file_extension(m_filename, ".fasta.gz")){
		return FASTA;
	}
	
	if( find_file_extension(m_filename, ".fastq") || find_file_extension(m_filename, ".fastq.gz") ){
		return FASTQ;
	}
	
	return UNKNOWN_SEQUENCE;
}

void SequenceIterator::load(const string &m_filename)
{
	// Clean up any open file pointers or memory buffers
	clear();
	
	file_type = get_file_type(m_filename);
	
	switch(file_type){
	
		case FASTA:
		case FASTQ:
		
			// Use zlib to read both compressed and uncompressed fasta files.
			fin = gzopen(m_filename.c_str(), "r");

			if(fin == NULL){

				cerr << "Error opening: " << m_filename << endl;
				throw __FILE__ ":SequenceIterator::SequenceIterator: Unable to open sequence file";
			}
			
			break;
		default:
			throw __FILE__ ":SequenceIterator::SequenceIterator: Unknown file type";
	}
	
	next();
}

void SequenceIterator::next()
{
	switch(file_type){
		case FASTA:
			next_fasta();
			break;
		case FASTQ:
			next_fastq();
			break;
		default:
			throw __FILE__ ":SequenceIterator::next: Unknown sequence file type";
			break;
	};
}

void SequenceIterator::next_fasta()
{
	if(fin == NULL){
		return;
	}
	
	const int buffer_len = 2048;
	char buffer[buffer_len];

	deque<char> seq_buffer;
	deque<char> info_buffer;
	
	while( gzgets(fin, buffer, buffer_len) ){
		
		if(strchr(buffer, '>') != NULL){
			
			info_buffer.clear();
			
			for(char* p = buffer;*p != '\0';++p){
				
				if( (*p != '\n') && (*p != '\r') ){
					info_buffer.push_back(*p);
				}
			}
			
			// Is this defline longer than the buffer?
			if(strpbrk(buffer, "\n\r") == NULL ){
				
				// Keep reading until we find an end of line symbol
				while( gzgets(fin, buffer, buffer_len) && !strpbrk(buffer, "\n\r") ){
					for(char* p = buffer;*p != '\0';++p){
						if( (*p != '\n') && (*p != '\r') ){
							info_buffer.push_back(*p);
						}
					}
				}
			}
			
			// Remove the fasta defline delimeter
			while(!info_buffer.empty() && ( isspace(info_buffer[0]) || (info_buffer[0] == '>') ) ){
				info_buffer.pop_front();
			}
			
			if( !seq_buffer.empty() ){
	
				swap(curr_defline, next_defline);
				
				next_defline.assign( info_buffer.begin(), info_buffer.end() );
				
				seq.assign( seq_buffer.begin(), seq_buffer.end() );
				return;
			}
			else{
				next_defline.assign( info_buffer.begin(), info_buffer.end() );
			}			
		}
		else{

			for(char* p = buffer;*p != '\0';++p){
				
				if( !isspace(*p) ){
				
					// Make sure all bases are upper-case
					seq_buffer.push_back( toupper(*p) );
				}
			}
		}
	}
	
	if( !seq_buffer.empty() ){
		
		swap(curr_defline, next_defline);
		
		seq.assign( seq_buffer.begin(), seq_buffer.end() );
		return;
	}
	
	gzclose(fin);
	fin = NULL;
}

void SequenceIterator::next_fastq()
{
	if(fin == NULL){
		return;
	}
	
	const int buffer_len = 2048;
	char buffer[buffer_len];
		
	deque<char> seq_buffer;
	deque<char> info_buffer;
	
	// The FASTQ file format is as follows for each record 
	// (from the Wikipedia entry):
	// Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
	// Line 2 is the raw sequence letters.
	// Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
	// Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence.

	/////////////////////////////////////////////////////////////////////////////////
	// Read the header
	/////////////////////////////////////////////////////////////////////////////////	
	while(true){
	
		if(gzgets(fin, buffer, buffer_len) == NULL){
			
			gzclose(fin);
			fin = NULL;
			
			return;
		}
		
		for(char* p = buffer;*p != '\0';++p){
		
			if( (*p != '\n') && (*p != '\r') ){
				info_buffer.push_back(*p);
			}
		}
		
		// If we found an end of line symbol, we can stop reading
		if(strpbrk(buffer, "\n\r") != NULL){
			break;
		}
	}
	
	// Remove the fastq defline delimeter
	while(!info_buffer.empty() && ( isspace(info_buffer[0]) || (info_buffer[0] == '@') ) ){
		info_buffer.pop_front();
	}

	curr_defline.assign( info_buffer.begin(), info_buffer.end() );
	
	/////////////////////////////////////////////////////////////////////////////////
	// Read the sequence
	/////////////////////////////////////////////////////////////////////////////////
	while(true){
	
		if(gzgets(fin, buffer, buffer_len) == NULL){
			throw __FILE__ ":next_fastq: Unable to read sequence";
		}

		for(char* p = buffer;*p != '\0';++p){
				
			if( !isspace(*p) ){
			
				// Make sure all bases are upper-case
				seq_buffer.push_back( toupper(*p) );
			}
		}

		// If we found an end of line symbol, we can stop reading
		if(strpbrk(buffer, "\n\r") != NULL){
			break;
		}
	}
	
	// Skip the '+'. Note that we *don't* actually test to make sure we read a single '+' on a line
	// by itself.
	if(gzgets(fin, buffer, buffer_len) == NULL){
		throw __FILE__ ":next_fastq: Unable to read '+'";
	}
	
	if( strpbrk(buffer, "\n\r") == NULL ){
		throw __FILE__ ":next_fastq: Error reading '+' delimiter";
	}

	/////////////////////////////////////////////////////////////////////////////////
	// Read and discard the quality
	/////////////////////////////////////////////////////////////////////////////////	
	while(true){
	
		if(gzgets(fin, buffer, buffer_len) == NULL){
			throw __FILE__ ":next_fastq: Unable to read quality";
		}

		// If we found an end of line symbol, we can stop reading
		if(strpbrk(buffer, "\n\r") != NULL){
			break;
		}
	}
		
	if( !seq_buffer.empty() ){
		
		seq.assign( seq_buffer.begin(), seq_buffer.end() );
		return;
	}
	
	gzclose(fin);
	fin = NULL;
}
