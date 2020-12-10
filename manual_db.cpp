#include <iostream>
#include <fstream>
#include <sstream>

#include <stdlib.h>
#include <getopt.h>

#include "maestro.h"
#include "caldera.h"

using namespace std;

int main(int argc, char *argv[])
{		
	try{
		
		// Command line arguments
		// -d <input database file to read accessions from>
		// -s <status file to update>
		// --meta <metadata file to read>
		
		string database_filename;
		string status_filename;
		string metadata_filename;

		const char* options = "d:s:h?";
		int config_opt = 0;
		int long_index = 0;

		struct option long_opts[] = {
			{"meta", true, &config_opt, 1},
			{0,0,0,0} // Terminate options list
		};

		int opt_code;
		opterr = 0;

		bool print_usage = (argc == 1);
	
		while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

			switch( opt_code ){
				case 0:

					if(config_opt == 1){ // bits
							
						metadata_filename = optarg;
						break;
					}

					cerr << "Unknown flag!" << endl;
					break;
				case 'd':
					database_filename = optarg;
					break;
				case 's':
					status_filename = optarg;
					break;
				case 'h':
				case '?':
					print_usage = true;
					break;
				default:
					cerr << '\"' << (char)opt_code << "\" is not a valid option!" << endl;
					break;
			};
		}

		if(print_usage){

			cerr << "Usage:" << endl;
			cerr << "\t-d <input database file to read accessions from>" << endl;
			cerr << "\t-s <status file to update>" << endl;
			cerr << "\t--meta <metadata file to read>" << endl;
			return EXIT_SUCCESS;
		}

		if( database_filename.empty() ){

			cerr << "Please specify a database file to read" << endl;
			return EXIT_SUCCESS;
		}

		if( status_filename.empty() ){

			cerr << "Please specify a status file to update" << endl;
			return EXIT_SUCCESS;
		}

		if( metadata_filename.empty() ){

			cerr << "Please specify a metadata file to read" << endl;
			return EXIT_SUCCESS;
		}

		// Make sure we can open the database file for reading
		ifstream fdb(database_filename.c_str(), ios::binary);

		if(!fdb){

			cerr << "Unable to open the database file " << database_filename << " for reading accessions" << endl;
			return EXIT_FAILURE;
		}

		// Read the binary metadata file and extract:
		//	- The list of SRA accessions
		//	- The starting location of the each SRA accession metadata record in the metadata file
		// Note that the resulting list of (accession, location) will be sorted by accession for
		// fast lookup by accession
		vector< pair<SraAccession, size_t /*location*/> > accession_loc;

		if( parse_accession_loc(accession_loc, metadata_filename, true /*verbose*/) == false ){

			cerr << "Error reading metadata information from " << metadata_filename << endl;
			return EXIT_FAILURE;
		}

		const size_t num_sra = accession_loc.size();

		if(num_sra == 0){

			cerr << "Did not read any SRA accessions from the input metadata file" << endl;
			return EXIT_FAILURE;
		}

		unsigned char *status = new unsigned char [num_sra];

		if(!status){
			throw __FILE__ ":main: Unable to allocate status buffer";
		}

		memset(status, STATUS_INIT, num_sra);

		size_t database_index = 1; // <-- Start counting database files from 1

		// Read the status file
		if( restore_status(status_filename, status, num_sra, database_index, false /*create missing*/) == false){

			if(status != NULL){

				delete [] status;
				status = NULL;
			}	
		}

		// Read the list of accessions from the database file
		DBFileHeader header;
			
		binary_read(fdb, header);
		
		if(!fdb){
		
			cerr << "Unable to read database header" << endl;
			return EXIT_FAILURE;
		}
		
		cerr << "Header information for " << database_filename << endl;
		cerr << "\tmagic = " << header.magic << endl;
		cerr << "\tversion = " << header.version << endl;
		cerr << "\tcrc32 = " << std::hex << header.crc32 << std::dec << endl;
		cerr << "\tkmer_len = " << header.kmer_len << endl;
		cerr << "\tnum_hash = " << header.num_hash << endl;
		
		const size_t filter_len = header.filter_len();
		
		cerr << "\tfilter_len = " << filter_len << endl;
		
		cerr << "\tlog_2_filter_len = " << header.log_2_filter_len << endl;
		cerr << "\tnum_filter = " << header.num_filter << endl;
		
		switch(header.hash_func){
			case MURMUR_HASH_32:
				cerr << "\thash_func = Murmur32" << endl;
				break;
			case UNKNOWN_HASH:
				cerr << "\thash_func = Unknown" << endl;
				break;
			default:
				cerr << "\thash_func = Invalid" << endl;
				break;
		};
		
		// Jump to the start of annotation information (after the array of FilerInfo
		// locations, which we do not need since we'll be reading the FileInfo structures
		// sequentially)
		fdb.seekg( header.info_start + header.num_filter*sizeof(unsigned long int) );

		size_t num_success = 0;
		size_t num_fail = 0;

		cerr << "Updating: ";
		string info_buffer;

		for(unsigned int i = 0;i < header.num_filter;++i){

			FilterInfo info;

			binary_read(fdb, info);

			if(!fdb){

				cerr << "Error reading FilterInfo for filter " << i << endl;
				throw __FILE__ ":main: Error reading FilterInfo from database";
			}

			if(info.run_accession == INVALID_ACCESSION){

				cerr << "Warning: FilterInfo " << i << " has an invalid run accession!" << endl;
				++num_fail;
				continue;
			}

			for(string::const_iterator j = info_buffer.begin();j != info_buffer.end();++j){
				cerr << '\b';
			}

			for(string::const_iterator j = info_buffer.begin();j != info_buffer.end();++j){
				cerr << ' ';
			}

			for(string::const_iterator j = info_buffer.begin();j != info_buffer.end();++j){
				cerr << '\b';
			}

			info_buffer = accession_to_str(info.run_accession);

			cerr << info_buffer;

			try{
				const size_t index = get_accession_index(info.run_accession, accession_loc);

				status[index] = STATUS_DATABASE_SUCCESS;

				++num_success;
			}
			catch(...){

				cerr << "Unable to find a valid status file index for SRA accession " 
					<< accession_to_str(info.run_accession) << endl;
				
				++num_fail;
				break;
			}
		}

		cerr << endl;

		fdb.close();

		cerr << "Successfully updated the status of " << num_success << " out of " 
			<< header.num_filter << " accessions" << endl;


		if(num_success != header.num_filter){
			cerr << "Failed to update all accession -- the status file will not be modified" << endl;
		}
		else{

			if( !write_status(status_filename, status, num_sra, database_index) ){
				cerr << "** Warning! Unable to commit status information to disk!" << endl;
			}
			else{
				cerr << "Successfully updated status file with " << num_success 
					<< " accessions from " << database_filename << endl;
			}
		}

		if(status != NULL){

			delete [] status;
			status = NULL;
		}
	}
	catch(const char *error){

        cerr << "Caught the error: " << error << endl;
        return EXIT_FAILURE;
	}
	catch(...){

		cerr << "Caught an unhandled error" << endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
