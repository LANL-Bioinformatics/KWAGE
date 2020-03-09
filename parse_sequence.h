#ifndef __PARSE_SEQUENCE
#define __PARSE_SEQUENCE

#include <string>
#include <zlib.h> // For compressed text-based files

// The allowed sequence file types
typedef enum {
	FASTA,
	FASTQ,
	UNKNOWN_SEQUENCE
} FileType;
	
// A class for iterating through a sequence file, one sequence at a time
class SequenceIterator
{
	private:
		gzFile fin;
		FileType file_type;
		std::string seq;
		std::string curr_defline;
		std::string next_defline;
		
		void next();
		void next_fasta();
		void next_fastq();
	public:
		SequenceIterator()
		{
			file_type = FASTA;
			fin = NULL;
		};
		
		SequenceIterator(const std::string &m_filename)
		{
			fin = NULL;
			
			load(m_filename);
		};
		
		~SequenceIterator()
		{
			clear();
		};
		
		void load(const std::string &m_filename);
		
		inline void clear()
		{
			if(fin != NULL){
			
				gzclose(fin);
				fin = NULL;
			}
		};
		
		operator bool() const
		{
			switch(file_type){
				case FASTA:
				case FASTQ:
					return (fin != NULL);
				default:
					throw __FILE__ ":SequenceIterator:: bool(): Unknown sequence file type";
					break;
			};
			
			// We should never get here
			return false;
		};
		
		// Read the next sequence
		void operator++()
		{
			next();
		};
		
		const std::string& get_seq() const
		{
			return seq;
		};
		
		const std::string& get_info() const
		{
			return curr_defline;
		};
};

#endif // __PARSE_SEQUENCE
