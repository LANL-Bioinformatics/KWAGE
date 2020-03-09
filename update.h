#ifndef __UPDATE
#define __UPDATE

#include <string>
#include <sstream>
#include <fstream>

class UpdateInfo : public std::stringstream {
	
	private:
		bool valid;
		size_t buffer_size;
		std::ofstream flog;
		
	public:
	
		UpdateInfo(bool m_valid, const std::string &m_log_file) :
			valid(m_valid), buffer_size(0)
		{
			log(m_log_file);
		};
		
		UpdateInfo(bool m_valid) :
			valid(m_valid), buffer_size(0)
		{
		};
		
		void log(const std::string &m_log_file)
		{
			if(!valid){
				return;
			}
			
			// Only attempt to open a log if the user has provided
			// a non-empty string
			if( !m_log_file.empty() ){
			
				flog.open( m_log_file.c_str() );

				if(!flog){
					throw __FILE__ ":Log: Unable to open the log file";
				}
			}
		};
		
		void flush();
		void close();
};

#endif // __UPDATE
