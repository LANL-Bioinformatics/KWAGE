// Download SRA sequence and metadata
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <deque>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#include <stdlib.h>
#include <getopt.h>
#include <signal.h>
#include <math.h>

#include "options.h"
#include "bloom.h"
#include "hash.h"
#include "word.h"
#include "parse_tar.h"
#include "parse_sequence.h"
#include "file_util.h"
#include "keys.h"
#include "bigsi++.h" // For DOWNLOAD_VERSION

using namespace std;

typedef enum{
	RUN_XML,
	EXPERIMENT_XML,
	SAMPLE_XML,
	STUDY_XML,
	//SUBMISSION_XML,
	//ANALYSIS_XML,
	//SRA_RUN_MEMBERS,
	SRA_ACCESSIONS,
	UNKNOWN
} SRAFileType;

#define		KB		size_t(1024)
#define		MB		(size_t(1024)*KB)
#define		GB		(size_t(1024)*MB)
#define		SEC_PER_DAY	(60*60*24)

bool match_prefix(const string &m_query, const string &m_subject);
string extract_prefix(const string &m_subject);
SRAFileType sra_file_type(const string &m_filename);
string parse_accession(const string &m_buffer);
string parse_xml(const string &m_key, const string &m_buffer);
void parse_sra_metadata(MAP<string, FilterInfo> &m_db, const string &m_metadata_file);
unordered_set<string> parse_sra_accessions(const string &m_logfile);
bool has_aligned_reads(const string &m_dir);

extern const char* allowed_filter_extention;

int main(int argc, char* argv[])
{
	DownloadOptions opt(argc, argv);

	if(opt.quit){
		return EXIT_SUCCESS;
	}

	time_t profile = time(NULL);

	// Before we start logging, we need to parse the log file (if it exists), to
	// extract any SRA accession that have already been downloaded
	const unordered_set<string> existing_sra_accessions = parse_sra_accessions(opt.log_file);
	
	cerr << "Found " << existing_sra_accessions.size() 
		<< " previously downloaded SRA records to skip" << endl;

	ofstream log(opt.log_file.c_str(), ios::app);

	log << "Started logging: " << ctime(&profile);

	// Save the command line arguments in the log file
	log << "Command-line (v. " << DOWNLOAD_VERSION << "):";
	
	for(int i = 0;i < argc;++i){
		log << ' ' << argv[i];
	}
	
	log << endl;
	
	// If we will be downloading, make sure that the scratch and download directories exist
	if(!opt.list_only){

		if( !make_dir(opt.download_dir) ){
			
			cerr << "Unable to create download directory: " << opt.download_dir << endl;
			throw __FILE__ ":main: Unable to create download directory";
		}
		
		if( !make_dir(opt.bloom_dir) ){
			
			cerr << "Unable to create Bloom filter directory: " << opt.bloom_dir << endl;
			throw __FILE__ ":main: Unable to create Bloom filter directory";
		}
	}

	// Step 1: Parse the SRA metadata to extract all of the SRA runs and selected metadata
	MAP<string, FilterInfo> db;

	cerr << "Reading NCBI metadata from " << opt.metadata_file << endl;
	log << "Reading NCBI metadata from " << opt.metadata_file << endl;

	time_t profile_parse_xml = time(NULL);

	parse_sra_metadata(db, opt.metadata_file);

	profile_parse_xml = time(NULL) - profile_parse_xml;
		
	const size_t num_sra_run = db.size();

	cerr << "Found metadata for " << num_sra_run << " SRA runs in " 
		<< profile_parse_xml << " sec" << endl;
	
	log << "Found metadata for " << num_sra_run << " SRA runs in " 
		<< profile_parse_xml << " sec" << endl;

	if(opt.list_only){

		for(MAP<string, FilterInfo>::const_iterator i = db.begin();i != db.end();++i){
			
			cout << i->first << endl;
			cout << "\texperiment_accession : " << i->second.experiment_accession << endl;
			cout << "\texperiment_title : " << i->second.experiment_title << endl;
			cout << "\texperiment_design_description : " << i->second.experiment_design_description << endl;
			cout << "\texperiment_library_name : " << i->second.experiment_library_name << endl;
			cout << "\texperiment_library_strategy : " << i->second.experiment_library_strategy << endl;
			cout << "\texperiment_library_source : " << i->second.experiment_library_source << endl;
			cout << "\texperiment_library_selection : " << i->second.experiment_library_selection << endl;
			cout << "\texperiment_instrument_model : " << i->second.experiment_instrument_model << endl;
			cout << "\tsample_accession : " << i->second.sample_accession << endl;
			cout << "\tsample_taxa : " << i->second.sample_taxa << endl;

			if( !i->second.sample_attributes.empty() ){

				cout << "\tsample_attributes" << endl;

				for(MULTIMAP<string, string>::const_iterator j = i->second.sample_attributes.begin();
					j != i->second.sample_attributes.end();++j){
					
					cout << "\t\t" << j->first << " : " << j->second << endl;
				}
			}

			cout << "\tstudy_accession : " << i->second.study_accession << endl;
			cout << "\tstudy_title : " << i->second.study_title << endl;
			cout << "\tstudy_abstract : " << i->second.study_abstract << endl;
		}

		return EXIT_SUCCESS;
	}

	cerr << "Extracting database keys" << endl;
	
	const vector<string> accessions = keys(db);
	
	const size_t num_accessions = accessions.size();
	
	cerr << "Found " << num_accessions 
		<< " SRA accession to download.\nDownloading sequence data using "
		<< opt.num_download_threads << " threads" << endl;
	
	time_t profile_download = time(NULL);
	
	size_t curr_run = 0;

	// The current directory is needed to generate the scripts that will
	// be submitted to the scheduler for Bloom filter creation.
	char* working_dir = getcwd(NULL, 0);

	if(working_dir == NULL){
		throw __FILE__ ":main_leader: Unable to obtain working directory";
	}
	
	// Reduce clutter for informative messages
	cerr << setprecision(3);
	log << setprecision(3);

	size_t total_sra_download_size = 0;
	
	#pragma omp parallel for num_threads(opt.num_download_threads)
	for(size_t i = 0;i < num_accessions;++i){
	//for(size_t i = 0;i < 2;++i){
	
		MAP<string, FilterInfo>::const_iterator db_iter = db.find(accessions[i]);
		
		if( db_iter == db.end() ){
			throw __FILE__ ":main: Unable to look up SRA accession";
		}

		const bool skip_download = 
			existing_sra_accessions.find(db_iter->first) != existing_sra_accessions.end();
		
		#pragma omp critical
		{
			++curr_run;

			cerr << "Downloading " << db_iter->first << " (" << curr_run << "): "
				<< (100.0*curr_run)/num_sra_run << '%' << endl;

			if(skip_download){
				cerr << "Skiping " << db_iter->first << " (already downloaded)" << endl;
			}
		}
		
		if(skip_download){
			continue;
		}
		
		// Change the working directory to the current run directory. This is needed
		// to ensure that any "helper" files that may be associated with SRA run
		// files (i.e. ".sra.vdbcache" and external reference sequences that are
		// required to decode aligned reads) are contained *within* the run directory.
		// This is will make it easier to clean up!
		const string run_dir = opt.download_dir + PATH_SEPARATOR + db_iter->first + "_run";
		
		if( !make_dir(run_dir) ){
			throw __FILE__ ":main_leader: Unable to create run directory";
		}

		// Note that the working directory is *process* dependent on Linux -- not thread
		// dependent! As a result, we cannot chdir() within OpenMP (or pthread) threads.
		//chdir( run_dir.c_str() );
		
		// Write the metadata to a file in the run_dir
		const string meta_name = run_dir + PATH_SEPARATOR + db_iter->first + ".info";
		
		ofstream fmeta(meta_name.c_str(), ios::binary);
		
		if(!fmeta){
			throw __FILE__ ":main: Unable to write metadata info file";
		}
		
		binary_write(fmeta, db_iter->second);
		
		fmeta.close();
		
		// For initial development only -- will need to put prefetch in the path
		// ** Prefetch will (by default) create a subdirectory with the same name
		// ** as the run accession (trying to change this with -O and/or -o
		// ** was not successfull)
		//
		// Note that the "prefetch" executable is actually a perl script that wraps
		// the executable shown below.
		//
		// As of Jan 27, 2020, NCBI has informed me that new versions of the SRA toolkit
		// will not support Aspera downloads. Hence the use of the "-t http" (even though
		// prefetch should figure this out on its own).
		//
		// The default max-size for prefetch is 20G; After increasing it to 100G, I have observed
		// that there are some SRA records much larger than 20G that yeild Bloom filter sizes that are
		// *greater* than the current maximum allowed Bloom filter size (currently 2GB).
		// 	55G	SRR10054704
		// 	54G	SRR10054705
		//	68G	SRR5738871
		//	69G	SRR5738872

		const size_t max_download_size = 30*GB;
		
		stringstream fetch_cmd;
		
		fetch_cmd << "$HOME/src/BIGSI/SRA/bin/prefetch-orig.2.10.0 --max-size " 
			<< max_download_size << " -t http -O " << run_dir << " " << db_iter->first;
		
		// Track the number of download attempts
		unsigned int num_attempts = 1;

		time_t fetch_profile = time(NULL);
		
		while(num_attempts <= opt.max_num_download_attempts){

			++num_attempts;

			system( fetch_cmd.str().c_str() );
		
			// If the SRA download was successfull, a run/run.sra file
			// will have been created in the current directory.
			const string sra_name = run_dir + PATH_SEPARATOR + 
				db_iter->first + PATH_SEPARATOR + 
				db_iter->first + ".sra";
			
			if( is_file(sra_name) ){
				break;
			}

			cerr << "\tRetrying download (attempt " << num_attempts << ")" << endl;
		}
		
		if(num_attempts > opt.max_num_download_attempts){

			cerr << "\tFailed to download " << db_iter->first << " after " << (num_attempts - 1) 
				<< " attempts" << endl;
			log << "Failed to download " << db_iter->first << " after " << (num_attempts - 1) 
				<< " attempts" << endl;
			
			// Remove any files or directories left over from this failed attempt
			// For example, attempting to download a controlled access SRA records
			// (with data from dbgap) will create a directory and metadata file, but
			// no SRA data.
			remove_all(run_dir);
			
			continue;
		}
		
		// For SRA records with associated alignments, the "prefetch" tool likes to create a subdirectory
		// in the *current working directory* to store the reference sequences. These references *should* be
		// placed in: run_dir + PATH_SEPARATOR + db_iter->first.
		// I suspect this is a bug (not a feature!) in prefetch. However, until it is fixed, we need to manually test
		// for the existance of a subdirectory with the current accession name, and, if it exists, move the contents
		// to the correct location!
		if( is_dir(db_iter->first) ){
			
			cerr << "\tMoving SRA associated files to the correct location" << endl;
			
			// Please note that the current implementation of move_files uses the "rename()" library function
			// to move files -- this will only work if both the source and destination directories are on
			// the same file system!
			const bool ret = move_files(db_iter->first, // source directory
				run_dir + PATH_SEPARATOR + db_iter->first, // destination directory
				true); // Remove source directory
			
			#pragma omp critical
			if(!ret){
			
				cerr << "\tError moving SRA associated files" << endl;
				log << "Error moving SRA associated files for " << db_iter->first << endl;
			}
		} 
		
		fetch_profile = time(NULL) - fetch_profile;
		
		const string sra_file_name = run_dir + PATH_SEPARATOR + db_iter->first + 
			PATH_SEPARATOR + db_iter->first + ".sra";
		const size_t sra_size = file_size(sra_file_name);
		
		#pragma omp critical
		{
			cerr << "Downloaded " << db_iter->first << " in " << fetch_profile << " sec ("
				<< double(sra_size)/GB << " GB)" << endl;

			////////////////////////////////////////////////////////////////////////////////
			// "Downloaded" Is a keyword in the log file that we will use when restarting an
			// SRA download to skip already downloaded records.
			log << "Downloaded " << db_iter->first  << endl;
			////////////////////////////////////////////////////////////////////////////////
		}
		
		// Check to see if the current SRA data directory contains a ".vdbcache" file.
		// If so, we assume that the SRA files contain aligned reads that are more 
		// computationally intensive to decompress. To handle this increased computational
		// burden, use additional threads to read and decompress the SRA data.
		const bool is_aligned_reads = has_aligned_reads(run_dir);
		
		// Use the cluster queueing to schedule the conversion of the SRA data into a Bloom filter
		// using a script that we will generate a run-time and pipe to qsub command
		stringstream bloom_cmd;
		
		// 1) Note that we need to escape the '!' to keep bash from try to interpret the command
		//    as an event from the history buffer
		// 2) The 'working_dir' string is the directory that we executed the sra_download program in. It
		//    is also used as the home directory for the cluster script to make sure that relative and
		//    absolute paths are processed correctly.
		// 3) The MPI parameters (--bind-to socket --oversubscribe) were hand-tuned for good
		//    performance on beagle.lanl.gov. Other (i.e. newer) machines will likely need additional
		//    tuning.
		// 4) The number of threads to use for reading the SRA data file (5 for "regular" reads and
		//    11 for aligned reads) is based on manual benchmarking on beagle.lanl.gov.
		//    If we use too many threads to read "regular" (i.e. unaligned) reads, we end up overwhelming
		//    the NFS file system (might work if we ever get a parallel file system).
		//	- The values of 5 (not 6) and 11 (not 12) are choosen to account for the extra thread
		//        spawned by the API for reading SRA files.
		// 5) Disabled the infinband network (using "--mca btl ^openib") until I can figure out why 
		//    it is not working with OpenMPI v 3 (I may just need to update beagle!!).
		
		const string bloom_file_name = opt.bloom_dir + PATH_SEPARATOR 
			+ db_iter->first + allowed_filter_extention;
		
		// Need to make these user configurable
		const size_t ranks_per_node = 4; // For beagle.lanl.gov
		const size_t bytes_per_node = size_t(10)*GB; // For beagle.lanl.gov
		
		const size_t num_node = max(size_t(1),
			size_t( ceil(double(sra_size)/bytes_per_node) ) );
		
		// Maximum number of attempts to create a Bloom filter
		const size_t max_bloom_iter = 3;
		
		bloom_cmd << "printf \"#\\!/bin/sh\\n\\n"
			<< "cd " << working_dir << "\\n"
			<< "iter=1\\n"
			<< "while [ \\$iter -le " << max_bloom_iter << " ] \\n"
			<< "do\\n"
			<< "mpirun -np " << ranks_per_node*num_node << " --bind-to socket --oversubscribe --mca btl ^openib "
			<< "\\$HOME/src/BIGSI/bigsi++1/bloomer -v -i " << run_dir
			<< " -o " << bloom_file_name
			<< " -k " << opt.kmer_len << " -p " << opt.false_positive_probability 
			<< " --len.min " << opt.min_log_2_filter_len
			<< " --len.max " << opt.max_log_2_filter_len
			<< " --slice " << (is_aligned_reads ? 11 : 5) << "\\n"
			<< "STATUS=\\$?\\n"
			<< "if [ \\$STATUS == 0 ] \\n"
			<< "then\\n"
			<< "iter=" << max_bloom_iter + 1 << "\\n" // Success! Break out of the loop!
			<< "else\\n"
			<< "iter=\\$(( \\$iter + 1 ))\\n" // Failure! Increment the loop counter
			<< "fi\\n"
			<< "done\\n"
			<< "\" | qsub -lnodes=" << num_node
			//<< " -e \\$HOME/src/BIGSI/bigsi++1/debug_status -o \\$HOME/src/BIGSI/bigsi++1/debug_status";
			<< " -e /dev/null -o /dev/null";
		
		//cerr << bloom_cmd.str() << endl;
		
		system( bloom_cmd.str().c_str() );
		
		#pragma omp atomic
		total_sra_download_size += sra_size;
		
		cerr << "Computing Bloom filter using " << ranks_per_node 
			<< " MPI ranks on " << num_node << " node" 
			<< ( (num_node == 1) ? "" : "s") << endl;
		cerr << "Current download rate is " 
			<< (double(total_sra_download_size)/GB)/(double(time(NULL) - profile_download)/SEC_PER_DAY)
			<< " GB/day" << endl;
			
		// Optionally sleep to (a) prevent downloading faster than we can convert
		// SRA files into Bloom filters and (b) avoid making the good folks at the NCBI
		// annoyed with us ...
		if(opt.sleep_interval > 0){
			sleep(opt.sleep_interval);
		}
		
		if(opt.max_backlog > 0){
			
			const unsigned int backlog_sleep = 30;
			size_t sleep_iter = 0;
			
			while( count_subdirectories(opt.download_dir) >= opt.max_backlog ){
				
				#pragma omp critical
				if(sleep_iter == 0){
					
					// Only log the first time
					cerr << "Maximum download backlog exceeded; sleeping" << endl;
					log << "Maximum download backlog exceeded; sleeping" << endl;
				}
				
				sleep(backlog_sleep);
				
				++sleep_iter;
			}
			
			#pragma omp critical
			if(sleep_iter > 0){
				log << "Slept for " << backlog_sleep*sleep_iter << " sec" << endl;
			}
		}		
	}
	
	if(working_dir != NULL){

		free(working_dir);
		working_dir = NULL;
	}

	profile_download = time(NULL) - profile_download;
	
	profile = time(NULL) - profile;
	
	cerr << "Finished download in " << profile_download << " sec" << endl;
	cerr << "Total run time was " << profile << " sec" << endl;
	cerr << "Total SRA download size was " << double(total_sra_download_size)/GB << " GB" << endl;
	
	log << "Finished download in " << profile_download << " sec" << endl;
	log << "Total run time was " << profile << " sec" << endl;
	log << "Total SRA download size was " << double(total_sra_download_size)/GB << " GB" << endl;

	profile = time(NULL);

	log << "Stopped logging: " << ctime(&profile);

    return EXIT_SUCCESS;
}

void parse_sra_metadata(MAP<string, FilterInfo> &m_db, const string &m_metadata_file)
{
	TarIterator iter(m_metadata_file);

	string curr_prefix = "/";
	string run_accession;
	string experiment_accession;
	string sample_accession;
	string sample_attribute_tag;
	string sample_attribute_value;
	bool sample_attribute = false;
	
	string study_accession;
	string study_title;
	string study_abstract;

	// Map accessions between different record types
	MAP<string, string> run_to_experiment;
	MAP<string, string> experiment_to_sample;
	MAP<string, string> experiment_to_title;

	MAP<string, string> experiment_to_design_description;
	MAP<string, string> experiment_to_library_name;
	MAP<string, string> experiment_to_library_strategy;
	MAP<string, string> experiment_to_library_source;
	MAP<string, string> experiment_to_library_selection;
	MAP<string, string> experiment_to_instrument_model;
	
	MAP<string, bool> experiment_is_controlled;
	MAP<string, string> sample_to_taxa;
	MULTIMAP< string, pair<string, string> > sample_to_attributes;
	
	deque<string> invalid_accessions;
	
	size_t line_number = 0;
	
	while(iter){
		
		++line_number;
		
		// DEBUG
		//if((*iter).find("controlled_access") != string::npos ){
		//	cerr << iter.filename() << " --> " << sra_file_type( iter.filename() ) << endl;
		//	cerr << '[' << iter.filename() << "] @ " << line_number << " : " << *iter << endl;
		//}
		
		//++iter;
		//continue;
		
		// Do we need to parse this file?
		bool valid_xml_file = true;
		
		switch( sra_file_type( iter.filename() ) ){
			case RUN_XML:
				
				// Append a space to the search key "RUN " to avoid picking up
				// "RUN_SET".
				if( (*iter).find("<RUN ") != string::npos ){
					run_accession = parse_accession(*iter);
				}
				
				if( (*iter).find("<EXPERIMENT_REF") != string::npos ){
				
					experiment_accession = parse_accession(*iter);
					
					if( run_accession.empty() ){
						throw __FILE__ ":main: Orphaned experiment accession";
					}
					
					run_to_experiment[run_accession] = experiment_accession;
				}
				
				break;
			case EXPERIMENT_XML:
				
				// Append a space to the search key "EXPERIMENT " to avoid picking up
				// "EXPERIMENT_SET".
				if( (*iter).find("<EXPERIMENT ") != string::npos ){
					experiment_accession = parse_accession(*iter);
				}
				
				// The keyword "<SAMPLE_DESCRIPTOR>" appears in SRA303418/SRA303418.experiment.xml.
				// To avoid examples like this, append a space to the end of the search string.
				if( (*iter).find("<SAMPLE_DESCRIPTOR ") != string::npos ){
				
					sample_accession = parse_accession(*iter);
					
					if( experiment_accession.empty() ){
						throw __FILE__ ":main: Orphaned sample accession";
					}
					
					experiment_to_sample[experiment_accession] = sample_accession;
				}
				
				if( (*iter).find("<TITLE>") != string::npos){
					
					const string title = parse_xml("TITLE", *iter);
					
					if( experiment_accession.empty() ){
						
						cerr << '[' << iter.filename() << "] @ line " << line_number << endl;
						throw __FILE__ ":main: Orphaned experiment title";
					}
					
					experiment_to_title[experiment_accession] = title;
				}
				
				if( (*iter).find("<DESIGN_DESCRIPTION>") != string::npos){
					
					const string design_description = parse_xml("DESIGN_DESCRIPTION", *iter);
					
					if( experiment_accession.empty() ){
						
						cerr << '[' << iter.filename() << "] @ line " << line_number << endl;
						throw __FILE__ ":main: Orphaned experiment design description";
					}
					
					experiment_to_design_description[experiment_accession] = design_description;
				}

				if( (*iter).find("<LIBRARY_NAME>") != string::npos){
					
					const string library_name = parse_xml("LIBRARY_NAME", *iter);
					
					if( experiment_accession.empty() ){
						
						cerr << '[' << iter.filename() << "] @ line " << line_number << endl;
						throw __FILE__ ":main: Orphaned experiment library name";
					}
					
					experiment_to_library_name[experiment_accession] = library_name;
				}

				if( (*iter).find("<LIBRARY_STRATEGY>") != string::npos){
					
					const string library_strategy = parse_xml("LIBRARY_STRATEGY", *iter);
					
					if( experiment_accession.empty() ){
						
						cerr << '[' << iter.filename() << "] @ line " << line_number << endl;
						throw __FILE__ ":main: Orphaned experiment library strategy";
					}
					
					experiment_to_library_strategy[experiment_accession] = library_strategy;
				}

				if( (*iter).find("<LIBRARY_SOURCE>") != string::npos){
					
					const string library_source = parse_xml("LIBRARY_SOURCE", *iter);
					
					if( experiment_accession.empty() ){
						
						cerr << '[' << iter.filename() << "] @ line " << line_number << endl;
						throw __FILE__ ":main: Orphaned experiment library source";
					}
					
					experiment_to_library_source[experiment_accession] = library_source;
				}

				if( (*iter).find("<LIBRARY_SELECTION>") != string::npos){
					
					const string library_selection = parse_xml("LIBRARY_SELECTION", *iter);
					
					if( experiment_accession.empty() ){
						
						cerr << '[' << iter.filename() << "] @ line " << line_number << endl;
						throw __FILE__ ":main: Orphaned experiment library selection";
					}
					
					experiment_to_library_selection[experiment_accession] = library_selection;
				}

				if( (*iter).find("<INSTRUMENT_MODEL>") != string::npos){
					
					const string instrument_model = parse_xml("INSTRUMENT_MODEL", *iter);
					
					if( experiment_accession.empty() ){
						
						cerr << '[' << iter.filename() << "] @ line " << line_number << endl;
						throw __FILE__ ":main: Orphaned experiment instrument model";
					}
					
					experiment_to_instrument_model[experiment_accession] = instrument_model;
				}

				// Check for the presence of dbgap data, to which access is controlled.
				// See: https://www.ncbi.nlm.nih.gov/books/NBK56558/
				if( (*iter).find("<EXTERNAL_ID namespace=\"dbgap\">") != string::npos){
					
					if( experiment_accession.empty() ){
						
						cerr << '[' << iter.filename() << "] @ line " << line_number << endl;
						throw __FILE__ ":main: Orphaned experiment dbgap id";
					}
					
					experiment_is_controlled[experiment_accession] = true;
				}
				
				break;
			case SAMPLE_XML:
				
				if( (*iter).find("<SAMPLE ") != string::npos ){
					sample_accession = parse_accession(*iter);
				}
				
				if( (*iter).find("<SCIENTIFIC_NAME>") != string::npos){
					
					const string taxa = parse_xml("SCIENTIFIC_NAME", *iter);
					
					if( sample_accession.empty() ){
						
						cerr << '[' << iter.filename() << "] @ line " << line_number << endl;
						throw __FILE__ ":main: Orphaned sample scientific name";
					}
					
					sample_to_taxa[sample_accession] = taxa;
				}
				
				if( (*iter).find("<SAMPLE_ATTRIBUTE>") != string::npos){
					sample_attribute = true;
				}
				
				if( (*iter).find("</SAMPLE_ATTRIBUTE>") != string::npos){
					sample_attribute = false;
				}
				
				if( sample_attribute && ( (*iter).find("<TAG>") != string::npos ) ){
					sample_attribute_tag = parse_xml("TAG", *iter);
				}
				
				if( sample_attribute && ( (*iter).find("<VALUE>") != string::npos ) ){
				
					sample_attribute_value = parse_xml("VALUE", *iter);
					
					if( sample_attribute_tag.empty() && sample_accession.empty() ){
						throw __FILE__ ":main: Orphaned sample attribute value";
					}
					
					// Skip the BioSampleModel tagged values (these entries are *not* reported
					// on the SRA BioSample web-pages).
					if(sample_attribute_tag != "BioSampleModel"){
					
						sample_to_attributes.insert( make_pair(sample_accession,
							make_pair(sample_attribute_tag, sample_attribute_value) ) );
					}
				}
				
				break;
			case STUDY_XML:

				if( (*iter).find("<PRIMARY_ID>") != string::npos){
					
					study_accession = parse_xml("PRIMARY_ID", *iter);
					
					if( study_accession.empty() ){
						
						cerr << '[' << iter.filename() << "] @ line " << line_number << endl;
						throw __FILE__ ":main: Orphaned study accession";
					}
				}

				if( (*iter).find("<PRIMARY_ID>") != string::npos){
					
					study_accession = parse_xml("PRIMARY_ID", *iter);
					
					if( study_accession.empty() ){
						
						cerr << '[' << iter.filename() << "] @ line " << line_number << endl;
						throw __FILE__ ":main: Orphaned study accession";
					}
				}

				if( (*iter).find("<STUDY_TITLE>") != string::npos){
					
					study_title = parse_xml("STUDY_TITLE", *iter);
					
					if( study_title.empty() ){
						
						cerr << '[' << iter.filename() << "] @ line " << line_number << endl;
						throw __FILE__ ":main: Orphaned study title";
					}
				}

				if( (*iter).find("<STUDY_ABSTRACT>") != string::npos ){
					
					study_abstract = parse_xml("STUDY_ABSTRACT", *iter);
					
					if( study_abstract.empty() ){
						
						cerr << '[' << iter.filename() << "] @ line " << line_number << endl;
						throw __FILE__ ":main: Orphaned study abstract";
					}
				}

				break;
			case SRA_ACCESSIONS:
			
				// The SRA_Accessions file does not count as a "valid" file as it does not
				// contain XML data.
				valid_xml_file = false;
				
				// While we can identify *some* invalid SRA experiments by looking for external ids
				// belonging to dbgap, this approach does not work for *all* SRA records. To catch
				// additional, non-public records we must also search the the SRA Accessions file
				
				//cerr << "In SRA Accession" << endl;
				
				if( ( (*iter).find("suppressed") != string::npos ) || 
				    ( (*iter).find("controlled_access") != string::npos ) ){
				    					
					// The accession (which could be run, experiment, sample, ...) is in the first
					// tab-delimited column
					const size_t loc = (*iter).find('\t');
				    
					if(loc != string::npos){
						
						const string accession = (*iter).substr(0, loc);
						
						invalid_accessions.push_back(accession);
					}
				}
				
				break;
			case UNKNOWN:
				
				valid_xml_file = false;
				break;
		};

		if(valid_xml_file == false){
			
			++iter;
			continue;
		}
		
		if( !match_prefix( curr_prefix, iter.filename() ) ){

			if( !run_accession.empty() ){
				
				// Output the required information for each run
				for(MAP<string, string>::const_iterator r = run_to_experiment.begin();
					r != run_to_experiment.end();++r){
					
					// Skip any controlled data
					const MAP<string, bool>::const_iterator is_controlled = experiment_is_controlled.find(r->second);
					
					if( ( is_controlled != experiment_is_controlled.end() ) && (is_controlled->second == true) ){
					
						// Don't write this to the terminal, since it will overwhelm the screen!
						//cerr << "Skipping " << r->first << " due to dbGaP access control rules" << endl;
						continue;
					}

					////////////////////////////////////////////////////////////////////////////////////
					// Create the metadata database entry
					FilterInfo &ref = m_db[r->first];
					
					if( !ref.experiment_accession.empty() && (ref.experiment_accession != r->second) ){
					
						cerr << "Warning! Experiment accession mismatch!" << endl;
						cerr << "\trun = " << r->first << endl;
						cerr << "\tinit experiment = " << ref.experiment_accession << endl;
						cerr << "\tnew experiment = " << r->second << endl;
					}
					
					ref.run_accession = r->first;
					ref.experiment_accession = r->second;
					
					const MAP<string, string>::const_iterator e2title = experiment_to_title.find(r->second);

					if( e2title != experiment_to_title.end() ){
						ref.experiment_title = e2title->second;
					}
					
					const MAP<string, string>::const_iterator e2design = experiment_to_design_description.find(r->second);

					if( e2design != experiment_to_design_description.end() ){
						ref.experiment_design_description = e2design->second;
					}

					const MAP<string, string>::const_iterator e2library_name = experiment_to_library_name.find(r->second);

					if( e2library_name != experiment_to_library_name.end() ){
						ref.experiment_library_name = e2library_name->second;
					}

					const MAP<string, string>::const_iterator e2library_strategy = experiment_to_library_strategy.find(r->second);

					if( e2library_strategy != experiment_to_library_strategy.end() ){
						ref.experiment_library_strategy = e2library_strategy->second;
					}

					const MAP<string, string>::const_iterator e2library_source = experiment_to_library_source.find(r->second);

					if( e2library_source != experiment_to_library_source.end() ){
						ref.experiment_library_source = e2library_source->second;
					}

					const MAP<string, string>::const_iterator e2library_selection = experiment_to_library_selection.find(r->second);

					if( e2library_selection != experiment_to_library_selection.end() ){
						ref.experiment_library_selection = e2library_selection->second;
					}

					const MAP<string, string>::const_iterator e2instrument_model = experiment_to_instrument_model.find(r->second);

					if( e2instrument_model != experiment_to_instrument_model.end() ){
						ref.experiment_instrument_model = e2instrument_model->second;
					}

					const MAP<string, string>::const_iterator x2s = experiment_to_sample.find(r->second);
					
					if( x2s == experiment_to_sample.end() ){
						//throw __FILE__ ":main: Unable to find sample for experiment";
						continue;
					}
					
					if( !ref.sample_accession.empty() && (ref.sample_accession != x2s->second) ){
					
						cerr << "Warning! Sample accession mismatch!" << endl;
						cerr << "\trun = " << r->first << endl;
						cerr << "\tinit sample = " << ref.sample_accession << endl;
						cerr << "\tnew sample = " << x2s->second << endl;
					}
					
					ref.sample_accession = x2s->second;
					
					const MAP<string, string>::const_iterator s2taxa = sample_to_taxa.find(x2s->second);

					if( s2taxa != sample_to_taxa.end() ){
					
						if( !ref.sample_taxa.empty() && (ref.sample_taxa != s2taxa->second) ){
							
							cerr << "Taxa name mismatch" << endl;
							cerr << "\trun = " << r->first << endl;
							cerr << "\tinit taxa = " << ref.sample_taxa << endl;
							cerr << "\tnew taxa = " << s2taxa->second << endl;
						}
						
						ref.sample_taxa = s2taxa->second;
					}
					
					typedef MULTIMAP< string, pair<string, string> >::const_iterator I;
					
					const pair<I, I> range = sample_to_attributes.equal_range(x2s->second);
					
					if(range.first != range.second){
					
						for(I i = range.first;i != range.second;++i){
							ref.sample_attributes.insert( make_pair(i->second.first, i->second.second) );
						}
					}

					ref.study_accession = study_accession;
					ref.study_title = study_title;
					ref.study_abstract = study_abstract;
				}
				
				// Clean up
				run_accession.clear();
				experiment_accession.clear();
				sample_accession.clear();
				sample_attribute_tag.clear();
				sample_attribute_value.clear();
				sample_attribute = false;
				
				run_to_experiment.clear();
				experiment_to_sample.clear();
				experiment_to_title.clear();
				experiment_to_design_description.clear();
				experiment_to_library_name.clear();
				experiment_to_library_strategy.clear();
				experiment_to_library_source.clear();
				experiment_to_library_selection.clear();
				experiment_to_instrument_model.clear();
				experiment_is_controlled.clear();
				sample_to_taxa.clear();
				sample_to_attributes.clear();

				study_accession.clear();
				study_title.clear();
				study_abstract.clear();
			}

			curr_prefix = extract_prefix( iter.filename() );
		}

		++iter;
		
		// DEBUG
		//if(m_db.size() >= 2000){
		//
		//	cerr << "DEBUG -- early parse termination!!" << endl;
		//	break;
		//}
	}

	cerr << "Found " << invalid_accessions.size() << " SRA run/experiment/sample accessions to exclude" << endl;
	
	size_t num_invalid = 0;
	
	for(deque<string>::const_iterator i = invalid_accessions.begin();i != invalid_accessions.end();++i){
		
		MAP<string, FilterInfo>::const_iterator iter = m_db.find(*i);
		
		if( iter != m_db.end() ){
		
			m_db.erase(iter);
			
			++num_invalid;
		}
	}
	
	cerr << "Removed " << num_invalid << " invalid/access controlled SRA runs" << endl;
}

bool match_prefix(const string &m_query, const string &m_subject)
{
	string::const_iterator q = m_query.begin();
	string::const_iterator s = m_subject.begin();

	while( ( q != m_query.end() ) && ( s != m_subject.end() ) ){
		
		if(*q != *s){
			return false;
		}

		++q;
		++s;
	}

	// If we got to the end of the prefix string, then we have a match
	return ( q == m_query.end() );
}

string extract_prefix(const string &m_subject)
{
	size_t loc = m_subject.find('/');

	if(loc == string::npos){
		return string();
	}

	return m_subject.substr(0, loc);
}

inline bool find_extension(const string &m_filename, const string &m_ext)
{
	const string::size_type loc = m_filename.rfind(m_ext);
	
	if(loc == string::npos){
		return false;
	}
	
	return ( ( loc + m_ext.size() ) == m_filename.size() );
}

// Match the file type by looking at the file name suffix
SRAFileType sra_file_type(const string &m_filename)
{
	
	if( find_extension(m_filename, ".run.xml") ){
		return RUN_XML;
	}
	
	if( find_extension(m_filename, ".experiment.xml") ){
		return EXPERIMENT_XML;
	}
	
	if( find_extension(m_filename, ".sample.xml") ){
		return SAMPLE_XML;
	}
	
	if( find_extension(m_filename, ".study.xml") ){
		return STUDY_XML;
	}
	
	if( find_extension(m_filename, "SRA_Accessions") ){
		return SRA_ACCESSIONS;
	}
	
	return UNKNOWN;
}

string parse_accession(const string &m_buffer)
{
	const size_t len = m_buffer.size();
	
	size_t loc = m_buffer.find("accession=\"");
	
	if(loc == string::npos){
		
		cerr << "m_buffer = " << m_buffer << endl;
		throw __FILE__ ":main: Unable to find key: accession=";
	}
	
	loc += 11; //strlen("accession=\"");
	
	for(size_t i = 1;(loc + i) < len;++i){
		
		if(m_buffer[loc + i] == '"'){
			return m_buffer.substr(loc, i);
		}
	}
	
	throw __FILE__ ":main: Unable to extract accession";
	
	return string();
}

string parse_xml(const string &m_key, const string &m_buffer)
{
	string query = "<" + m_key + ">";
	
	size_t begin = m_buffer.find(query);
	
	if(begin == string::npos){
		throw __FILE__ ":parse_xml: Unable to find begin key";
	}
	
	begin += query.size();
	
	query = "</" + m_key + ">";
	
	const size_t end = m_buffer.rfind(query);
	
	if(end == string::npos){
		throw __FILE__ ":parse_xml: Unable to find end key";
	}
	
	if(begin > end){
		throw __FILE__ ":parse_xml: XML delimeters are out of order";
	}
	
	return m_buffer.substr(begin, end - begin);
}

unordered_set<string> parse_sra_accessions(const string &m_logfile)
{
	unordered_set<string> ret;
	
	ifstream fin( m_logfile.c_str() );
	
	if(!fin){
		return ret;
	}
	
	const string keyword = "Downloaded ";

	string line;
	size_t line_number = 0;
	
	while( getline(fin, line) ){
		
		++line_number;
		
		if(line.find(keyword) == 0){
			
			line = line.substr( keyword.size(), line.size() - keyword.size() );
			
			if( line.empty() ){
				cerr << "** Warning ** Bad \"Download\" keyword at line " << line_number << endl;
			}
			else{
				ret.insert(line);
			}
		}
	}
	
	return ret;
}

bool has_aligned_reads(const string &m_dir)
{
	FindFiles ff(m_dir);
	
	while( ff.next() ){
		
		if( find_file_extension(ff.name(), ".vdbcache") ) {
			return true;
		}
		
	}
	
	return false;
}
