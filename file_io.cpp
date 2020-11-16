#include <algorithm>

#include "maestro.h"
#include "binary_io.h"
#include "file_util.h"

using namespace std;

size_t get_accession_index(const SraAccession &m_accession, 
	const vector< pair<SraAccession, size_t /*location*/> > &m_accession_loc)
{
	vector< pair<SraAccession, size_t /*location*/> >::const_iterator iter = 
		lower_bound( m_accession_loc.begin(), m_accession_loc.end(), m_accession, search_by_accession() );

	if( iter == m_accession_loc.end() ){
		throw __FILE__ ":get_accession_index: Unable to find accession";
	}

	return iter - m_accession_loc.begin();
}


bool parse_accession_loc(vector< pair<SraAccession, size_t /*location*/> > &m_accession_loc, const string &m_filename,
	bool m_verbose)
{
	ifstream fin(m_filename, ios::binary);

	if(!fin){

		cerr << "Unable to open metadata file (" << m_filename << ") for reading";
		return false;
	}

	size_t num_sra_meta = 0;

	binary_read(fin, num_sra_meta);

	if(!fin){

		cerr << "Unable to read the number of sra metadata records" << endl;
		return false;
	}

	if(m_verbose){
		cerr << "Reading " << num_sra_meta << " SRA records from the metadata file:" << endl;
	}

	m_accession_loc.resize(num_sra_meta);

	//string update;
	const size_t update_every = max(1UL, num_sra_meta/20);

	FilterInfo info;

	for(size_t i = 0;i < num_sra_meta;++i){

		// The start of the current annotation information for the i^th record
		const size_t loc = fin.tellg();

		try{
			binary_read(fin, info);
		}
		catch(...){

			cerr << "\nUnable to read SRA record " << i << " of " << num_sra_meta << endl;
			return false;
		}

		m_accession_loc[i] = make_pair(info.run_accession, loc);

		if( m_verbose && (i%update_every == 0) ){
		//	
		//	for(string::const_iterator j = update.begin();j != update.end();++j){
		//		cerr << '\b';
		//	}
		//
		//	for(string::const_iterator j = update.begin();j != update.end();++j){
		//		cerr << ' ';
		//	}
		//
		//	for(string::const_iterator j = update.begin();j != update.end();++j){
		//		cerr << '\b';
		//	}
		//
		//	stringstream ssout;
		//
		//	ssout << (100.0*i)/num_sra_meta << '%';
		//
		//	update = ssout.str();
		//
		//	cerr << update;
			// Round the % complete to the nearest integer
			cerr << ' ' << int( (100.0*i)/num_sra_meta ) << '%';
		}
	}

	// Sort by accession for fast look up
	sort( m_accession_loc.begin(), m_accession_loc.end() );

	if(m_verbose){
	//
	//	for(string::const_iterator j = update.begin();j != update.end();++j){
	//		cerr << '\b';
	//	}
	//
	//	for(string::const_iterator j = update.begin();j != update.end();++j){
	//		cerr << ' ';
	//	}
	//
	//	for(string::const_iterator j = update.begin();j != update.end();++j){
	//		cerr << '\b';
	//	}
	//	
		cerr <<  " done" << endl;
	}

	return true;
}

bool read_sra_repository(string &m_path)
{
	const string filename =  getenv("HOME") + string("/.ncbi/user-settings.mkfg");

	ifstream fin( filename.c_str() );

	if(!fin){
		
		cerr << "Unable to open SRA toolkit config file: " << filename << endl;
		return false;
	}

	const string repository_key = "/repository/user/main/public/root";

	string line;

	while( getline(fin, line) ){

		if(line.find(repository_key) != string::npos){

			size_t begin = repository_key.size();
			size_t end = line.size();

			// Skip leading white space
			while( (begin < end) && isspace(line[begin]) ){
				++begin;
			}

			// The next leading character must be a '='
			if( (begin >= end) || (line[begin] != '=') ){

				cerr << "Error parsing (missing '=') SRA toolkit config file: " << filename << endl;
				return false;
			}

			++begin;

			// Skip leading white space
			while( (begin < end) && isspace(line[begin]) ){
				++begin;
			}

			// The next leading character must be a '"'
			if( (begin >= end) || (line[begin] != '"') ){

				cerr << "Error parsing (missing '\"') SRA toolkit config file: " << filename << endl;
				return false;
			}

			++begin;

			// Skip trailing white space
			while( (begin < end) && isspace(line[end - 1]) ){
				--end;
			}

			// The next trailing character must be a '"'
			if( (begin >= end) || (line[end - 1] != '"') ){

				cerr << "Error parsing (missing '\"') SRA toolkit config file: " << filename << endl;
				return false;
			}

			--end;

			// The local SRA file repository directory has the following subdirectories:
			// 	files/	
			//	nannot/	
			//	refseq/	
			//	sra/	<-- Contains the downloaded .sra files
			//	wgs/
			m_path = line.substr(begin, end - begin) + PATH_SEPARATOR + string("sra");

			return true;
		}
	}

	cerr << "Unable to find '/repository/user/main/public/root' in SRA toolkit config file" << filename << endl;
	return false;
}

bool write_status(const string &m_filename, unsigned char* m_status, const size_t &m_num_sra, size_t &m_database_index)
{
	const string temp_filename = m_filename + ".temp";

	ofstream fout( temp_filename.c_str(), ios::binary );

	if(!fout){

		cerr << "Unable to create a status file" << endl;
		return false;
	}

	// Write the current (starting) database index
	fout.write( (char*)(&m_database_index), sizeof(size_t) );

	// Write the number of SRA accession elements
	fout.write( (char*)(&m_num_sra), sizeof(size_t) );
	fout.write( (char*)(m_status), m_num_sra);

	if(!fout){

		cerr << "Unable to write status file" << endl;
		return false;
	}

	// Move temp_filename to m_filename
	if(rename( temp_filename.c_str(), m_filename.c_str() ) != 0){

		cerr << "Error updating status file" << endl;
		return false;
	}

	return true;
}

bool restore_status(const string &m_filename, unsigned char* m_status, const size_t &m_num_sra, size_t &m_database_index,
	bool m_create_missing)
{
	ifstream fin( m_filename.c_str(), ios::binary );

	if(!fin){

		if(m_create_missing){

			// The status file does not exist, create it
			return write_status(m_filename, m_status, m_num_sra, m_database_index);
		}

		cerr << "Could not open status file " << m_filename << " for reading" << endl;

		return false;
	}

	// Read the database index
	fin.read( (char*)(&m_database_index), sizeof(size_t) );

	size_t len;

	fin.read( (char*)(&len), sizeof(size_t) );

	if(!fin){

		cerr << "Unable to read number of elements in existing status file" << endl;
		return false;
	}

	if(len != m_num_sra){

		cerr << "Error reading existing status file. Read " << len << ", but expected " << m_num_sra << endl;
		return false;
	}

	fin.read( (char*)m_status, m_num_sra );

	if(!fin){

		cerr << "Unable to read SRA status records from existing status file" << endl;
		return false;
	}

	return true;
}