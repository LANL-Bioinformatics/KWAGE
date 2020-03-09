#ifndef __FILE_UTIL
#define __FILE_UTIL

#include <string>
#include <deque>

#include <string.h> // strcmp
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <ctype.h>

#define	PATH_SEPARATOR	'/'

class FindFiles{

    private:
        
        // The target directories and/or files that remain to
        // be searched
        std::deque<std::string> targets;
        std::string curr_file;
        std::string curr_dir;

        DIR *dp;

    public:

        // Return true if curr_file has an entry, false otherwise
        bool next()
        {
            // Remove the last returned file (if any)
            curr_file.clear();

            if(dp == NULL){

                // Is there another directory to explore?
                if( targets.empty() ){
                    return false;
                }

                const std::string p = targets.front();
                
                targets.pop_front();

                // Is the input path a directory?
                struct stat path_info;

                if( stat(p.c_str(), &path_info) != 0){
                    throw __FILE__ ":FindFiles::next: Unable to stat entry";
                }

                if( S_ISREG(path_info.st_mode) ) { // A regular file

                    curr_file = p;
                    return true;
                }
                
                // Is this a directory
                if( S_ISDIR(path_info.st_mode) ){
                
                    curr_dir = p;

                    // Read the contents of the directory
                    dp = opendir( curr_dir.c_str() );

                    if(dp == NULL){
                        throw __FILE__ ":FindFiles::next: Unable to open directory for reading";
                    }

                    return next();
                }

                // If we get here, we have encountered a filesystem object that
                // is not a file or a directory!
                throw __FILE__ ":FindFiles::next: Unknown filesystem object";
            }

            struct dirent *dir = NULL;

            while( ( dir = readdir(dp) ) ){
                
                // Skip any removed files or diretories
                if(dir->d_ino == 0){
                    continue;
                }

                // Skip the special directories "." and ".."
                if( (strcmp(dir->d_name, ".") == 0) ||
                    (strcmp(dir->d_name, "..") == 0) ){
                    continue;
                }

                const std::string name = curr_dir + PATH_SEPARATOR + dir->d_name;
                        
                // Is this entry a directory or file name?
                struct stat dir_info;

                if( stat(name.c_str(), &dir_info) != 0){
                    throw __FILE__ ":FindFiles::next: Unable to stat entry (2)";
                }

                if( S_ISDIR(dir_info.st_mode) ){
                    
                    // Save the names of the subdirectories to search. Don't
                    // search them now to avoid the risk of opening too many
                    // directories for reading (since we do not know the system
                    // limit).
                    targets.push_back(name);
                }
                else if( S_ISREG(dir_info.st_mode) ){
                    
                        curr_file = name;
                        return true;
                }
            }

            closedir(dp);
            dp = NULL;
            curr_dir.clear();

            // We have exhausted the contents of the current directory.
            // Is there another directory to search?
            return next();
        };

        FindFiles() : dp(NULL)
        {
        };

        FindFiles(const std::string &m_path) : dp(NULL)
        {
            add(m_path);
        };

        template<class ITERATOR>
        FindFiles(const ITERATOR &m_begin, const ITERATOR &m_end) : dp(NULL)
        {
            add(m_begin, m_end);
        };

        ~FindFiles()
        {
            clear();
        };

        void add(const std::string &m_path)
        {
            targets.push_back(m_path);
        };

        template<class ITERATOR>
        void add(const ITERATOR &m_begin, const ITERATOR &m_end)
        {
            for(ITERATOR i = m_begin;i != m_end;++i){
                targets.push_back(*i);
            }
        };

        void clear()
        {
            targets.clear();
            curr_file.clear();
            curr_dir.clear();

            if(dp != NULL){

                closedir(dp);
                dp = NULL;
            }
        };

        inline const std::string& dir() const
        {
            return curr_dir;
        };

        inline const std::string& name() const
        {
            return curr_file;
        };
};

// This helper class is needed to sort directories by string length
// to ensure that child directories are deleted before parent directories
struct sort_by_length
{
	inline bool operator()(const std::string &m_a, const std::string &m_b)const
	{
		return m_a.size() < m_b.size();
	};
};

bool make_dir(const std::string &m_dirname);
bool is_dir(const std::string &m_dirname);
bool is_file(const std::string &m_filename);
size_t file_size(const std::string &m_filename);
bool match_extension(const std::string &m_input, const std::string &m_ext);
bool find_file_extension(const std::string &m_path, const char** m_ext);
bool find_file_extension(const std::string &m_path, const char* m_ext);
std::string strip_trailing_path_separator(const std::string &m_str);
size_t count_subdirectories(const std::string &m_dir);
std::pair<std::string /*dir*/, std::string /*file*/> 
	split_dir_and_file(const std::string &m_filename);
bool remove_all(const std::string &m_dir);
bool move_files(const std::string &m_src, const std::string &m_dst,
	const bool &m_rm_src);

#endif // __FILE_UTIL
