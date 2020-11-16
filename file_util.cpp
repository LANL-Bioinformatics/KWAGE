#include "file_util.h"
#include "ifind.h"

#include <algorithm>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>

using namespace std;

// Return true if the directory exists or if we can
// create it. Return false otherwise.
bool make_dir(const std::string &m_dirname)
{
	
	struct stat dir_info;

	// Read, write and execute for the owner, read-only for all others
	const mode_t default_mode = S_IRWXU | S_IRGRP | S_IROTH;
	
	if(stat(m_dirname.c_str(), &dir_info) != 0){
		
		// Try to create the directory
		if(mkdir(m_dirname.c_str(), default_mode) != 0){
			return false;
		}
		
		return true;
	}
	
	// The path exists, make sure it is a directory
	if( S_ISDIR(dir_info.st_mode) ){
		return true;
	}

	return false;
}

bool is_dir(const std::string &m_dirname)
{
	struct stat dir_info;

	if(stat(m_dirname.c_str(), &dir_info) != 0){
		return false;
	}
	
	// The path exists, make sure it is a directory
	return S_ISDIR(dir_info.st_mode);
}

bool is_file(const std::string &m_filename)
{
	struct stat file_info;

	if(stat(m_filename.c_str(), &file_info) != 0){
		return false;
	}
	
	// The path exists, make sure it is a file
	return S_ISREG(file_info.st_mode);
}

bool is_path(const string &m_name)
{
	struct stat path_info;

	if(stat(m_name.c_str(), &path_info) != 0){
		return false;
	}
	
	// The path exists, make sure it is a file
	return ( S_ISDIR(path_info.st_mode) || S_ISREG(path_info.st_mode) );
}

size_t file_size(const string &m_filename)
{
	struct stat file_info;

	if(stat(m_filename.c_str(), &file_info) != 0){
		return 0;
	}
	
	// The path exists, make sure it is a file
	if( S_ISREG(file_info.st_mode) ){
		return file_info.st_size;
	}

	return 0;
}

bool find_file_extension(const string &m_path, const char** m_ext)
{
	// For extension matching, we need case-insensitive string matching
	for(const char** ptr = m_ext;*ptr != NULL;++ptr){
		
		if( find_file_extension(m_path, *ptr) ){
			return true;
		}
	}
	
	return false;
}

bool find_file_extension(const string &m_path, const char* m_ext)
{
	// For extension matching, we need case-insensitive string matching
	size_t loc = ifind(m_path, m_ext);
		
	if(loc != string::npos){
			
		if( (strlen(m_ext) + loc) == m_path.size() ){
			return true;
		}
	}
	
	return false;
}

string strip_trailing_path_separator(const string &m_str)
{
	string ret(m_str);
	
	// Some older C++ compilers do not have string::pop_back()
	const size_t init_len = ret.size();
	size_t len = init_len;

	while( (len > 0) && ( ret[len - 1] == PATH_SEPARATOR ) ){
		--len;
	}

	if(len != init_len){
		ret = ret.substr(0, len);
	}
	
	return ret;
}

bool match_extension(const string &m_input, const string &m_ext)
{
	const size_t len = m_input.size();
	const size_t ext_len = m_ext.size();

	if(len < ext_len){
		return false;
	}

	return ( m_input.find(m_ext) == (len - ext_len) );
}

size_t count_subdirectories(const string &m_dir)
{
	size_t ret = 0;
	
	// Is the input path a directory?
	struct stat dir_info;

	if( stat(m_dir.c_str(), &dir_info) != 0){
		throw __FILE__ ":count_subdirectories: Unable to stat path";
	}

	// Is this a directory
	if( !S_ISDIR(dir_info.st_mode) ){
		throw __FILE__ ":count_subdirectories: Path is not a directory";
	}	

	// Read the contents of the directory
	DIR *dp = opendir( m_dir.c_str() );
	

	if(dp == NULL){
		throw __FILE__ ":count_subdirectories: Unable to open directory for reading";
	}

	struct dirent *d = NULL;

	while( ( d = readdir(dp) ) ){

		// Skip any removed files or diretories
		if(d->d_ino == 0){
			continue;
		}

		// Skip the special directories "." and ".."
		if( (strcmp(d->d_name, ".") == 0) ||
		    (strcmp(d->d_name, "..") == 0) ){
			continue;
		}

		
		const std::string name = m_dir + PATH_SEPARATOR + d->d_name;
		
		if( stat(name.c_str(), &dir_info) != 0){
                    throw __FILE__ ":count_subdirectories: Unable to stat entry (2)";
                }
		
		if( S_ISDIR(dir_info.st_mode) ){
			++ret;
		}
	}
	
	closedir(dp);
	
	return ret;
}

// If the filename has a path separator in it (i.e. "/"), return all characters
// *before* the last path separator as the directory and *after* the last path
// separator as the filename.
pair<string /*dir*/, string /*file*/>  split_dir_and_file(const string &m_filename)
{
	const size_t loc = m_filename.find_last_of(PATH_SEPARATOR);
	
	if(loc == string::npos){
		return pair<string, string>(string("."), m_filename);
	}
	
	return pair<string, string>( m_filename.substr(0, loc), 
		m_filename.substr(loc + 1, m_filename.size() - (loc + 1) ) );
}

string parent_dir(const string &m_filename)
{
	return split_dir_and_file(m_filename).first;
}

// /This/is/an/arbitrary/././example/path/foo.txt
//                                        ^^^^^^^
//                                         leaf
string leaf_path_name(const string &m_path)
{
	size_t end = m_path.size();

	if(end == 0){
		throw __FILE__ ":leaf_path_name: Empty path string (1)";
	}

	// Skip any terminal path separators
	while( (end > 0) && (m_path[end - 1] == PATH_SEPARATOR) ){
		--end;
	}

	if(end == 0){
		throw __FILE__ ":leaf_path_name: Empty path string (2)";
	}

	size_t begin = end - 1;

	while( (begin > 0) && (m_path[begin - 1] != PATH_SEPARATOR) ){
		--begin;
	}

	return m_path.substr(begin, end - begin);
}

// Return true if we were successfull in removing all files and directories, 
// false otherwise
bool remove_all(const string &m_dir)
{
	bool ret = true;
	
	FindFiles ff(m_dir);
		
	deque<string> del_dir;

	while( ff.next() ){

		if(unlink( ff.name().c_str() ) != 0){
			ret = false;
		}

		del_dir.push_back( ff.dir() );
	}

	del_dir.push_back(m_dir);

	// Make the list of directories to delete unique
	sort( del_dir.begin(), del_dir.end() );
	del_dir.erase( unique( del_dir.begin(), del_dir.end() ), del_dir.end() );

	sort( del_dir.begin(), del_dir.end(), sort_by_length() );

	for(deque<string>::const_reverse_iterator i = del_dir.rbegin();i != del_dir.rend();++i){

		if(rmdir( i->c_str() ) != 0){
			ret = false;
		}
	}
	
	return ret;
}

bool move_files(const string &m_src, const string &m_dst,
	const bool &m_rm_src)
{
	bool ret = true;
	
	FindFiles ff(m_src);

	while( ff.next() ){

		// Use rename to move files. This requires that both the source and destination
		// directories belong to the same file system. If this ends up being too restrictive,
		// we will need to use a system call (or write our own using read/write, etc)!
		const string dst_name = m_dst + PATH_SEPARATOR + split_dir_and_file( ff.name() ).second;
		
		if(rename( ff.name().c_str(), dst_name.c_str() ) != 0){
			ret = false;
		}
	}

	if(m_rm_src){
	
		if(rmdir( m_src.c_str() ) != 0){
			ret = false;
		}
	}
	
	return ret;

}

// Create any directories that do not already exist
bool create_path(const string &m_path)
{
	size_t loc = 0;

	while( ( loc = m_path.find(PATH_SEPARATOR, loc) ) != string::npos ){

		const string sub_path = m_path.substr(0, loc);

		if( !is_dir(sub_path) ){
			if( !make_dir(sub_path) ){
				return false;
			}
		}

		++loc;
	}

	return true;
}