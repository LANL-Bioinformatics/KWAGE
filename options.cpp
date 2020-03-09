#include <iostream> 
#include <unordered_set>
#include <unordered_map>
#include <deque>

#include <getopt.h>
#include <string.h>
#include <math.h>

#include "options.h"
#include "word.h"
#include "file_util.h"
#include "bigsi++.h"
#include "hash.h"

using namespace std;

// Enumerate the allowed fasta file extensions
const char* allowed_sequence_extentions [] = {
	".fna", ".fna.gz",
	".fasta", ".fasta.gz",
	".fa", ".fa.gz",
	".fastq", ".fastq.gz",
	".sra",
	NULL
};

const char* allowed_db_extentions [] = {
	".db",
	NULL
};

const char* allowed_filter_extention = ".bloom";

size_t str_to_uint32_t(const string &m_str);
bool is_fasta(const string &m_input);

void BuildOptions::load(int argc, char* argv[])
{	
	// Don't quit unless we need to
	quit = false;

	compress = false;
	subfilter_bits = DEFAULT_SUBFILTER_BITS;
	compression_threshold = DEFAULT_COMPRESSION_THRESHOLD;
	
	unordered_set<string> inputs;
	
	// Command line arguments
	// --write <output directory name>
	// --read <input Bloom filter directory> (can be repeated)
	// [--compress (compress Bloom filters in output files)]
	// [--compress.threshold <threshold> (Minimum fractional size reduction; [0,1])]
	// [--log <log filename>]

	const char* options = "h?";
	int config_opt = 0;
	int long_index = 0;

	struct option long_opts[] = {
		{"write", true, &config_opt, 2},
		{"read", true, &config_opt, 3},
		{"compress", false, &config_opt, 5},
		{"no-compress", false, &config_opt, 6},
		{"log", true, &config_opt, 7},
		{"compress.threshold", true, &config_opt, 8},
		{0,0,0,0} // Terminate options list
	};

	int opt_code;
	opterr = 0;

	bool print_usage = (argc == 1);

	while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

		switch( opt_code ){
			case 0:
				if(config_opt == 2){ // --write 
						
					output_dir = optarg;
					break;
				}
				
				if(config_opt == 3){ // --read 
						
					inputs.insert( strip_trailing_path_separator(optarg) );
					break;
				}
				
				if(config_opt == 5){ // --compress 
						
					compress = true;
					break;
				}
				
				if(config_opt == 6){ // --no-compress 
						
					compress = false;
					break;
				}
				
				if(config_opt == 7){ // --log 
						
					log_filename = optarg;
					break;
				}
				
				if(config_opt == 8){ // --compress.threshold 
						
					compression_threshold = atof(optarg);
					break;
				}
												
				cerr << "Unknown flag!" << endl;
				break;
			case '?':
			case 'h':
				print_usage = true;
				break;
			default:
				cerr << '\"' << (char)opt_code << "\" is not a valid option!" << endl;
				break;
		};
	}

	if(print_usage){

		quit = true;

		cerr << "Usage for BuildDB (v. " << BUILD_DB_VERSION << "):" << endl;
		cerr << "\t--write <output directory name>" << endl;
		cerr << "\t--read <input file directory> (can be repeated)" << endl;
		cerr << "\t[--compress (compress Bloom filters in output files)]" << endl;
		cerr << "\t[--compress.threshold <threshold> (Minimum fractional size reduction] (default is "
			<< DEFAULT_COMPRESSION_THRESHOLD << ")" << endl;
		cerr << "\t[--log <log filename>]" << endl;
		
		return;
	}

	if( inputs.empty() ){
	
		cerr << "Please specify one or more directories to search for Bloom filter files" << endl;
		quit = true;
		return;
	}
	
	if( (compression_threshold < 0.0) || (compression_threshold > 1.0) ){
		
		cerr << "Please specify: 0 < compression threshold <= 1.0" << endl;
		quit = true;
		return;
	}
	
	// Recursively search each input directory or file and extract all subdirectories that contain
	// Bloom filter files. 
	FindFiles ff( inputs.begin(), inputs.end() );

	while( ff.next() ){

		if( find_file_extension(ff.name(), allowed_filter_extention) ) {
			filter_input.push_back( ff.name() );
		}
	}

	if( filter_input.empty() ){

		cerr << "Unable to find an Bloom filter files" << endl;
			
		quit = true;
		return;
	}
}

SearchOptions::SearchOptions(int argc, char* argv[])
{	
	// Don't quit unless we need to
	quit = false;
	threshold = DEFAULT_SEARCH_THRESHOLD;
	output_format = DEFAULT_SEARCH_OUTPUT;
	
	// Command line arguments
	// [-o <output file name>]
	// [--o.csv (output CSV) | --o.json (output JSON)]
	// [-t <search threshold>] (default is 1)
	// -d <database search path> (can be repeated)
	// [-i <input sequence file>] (can be repeated)
	// [<DNA sequence>] (can be repeated)

	const char* options = "o:d:i:t:h?";
	int config_opt = 0;
	int long_index = 0;

	struct option long_opts[] = {
		{"o.csv", false, &config_opt, 1},
		{"o.json", false, &config_opt, 2},
		{0,0,0,0} // Terminate options list
	};

	int opt_code;
	opterr = 0;

	bool print_usage = (argc == 1);

	FindFiles ff;

	while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

		switch( opt_code ){
			case 0:
				if(config_opt == 1){ // o.csv 
						
					output_format = OUTPUT_CSV;
					break;
				}
				
				if(config_opt == 2){ // o.json 
						
					output_format = OUTPUT_JSON;
					break;
				}
				
				cerr << "Unknown flag!" << endl;
				break;
			case 'o':
				output_file = optarg;
				break;
			case 'i':
				query_files.push_back(optarg);
				break;
			case 'd':
				
				// Add this directory (or file) to the list of 
				// database files to search
				ff.add(optarg);
				break;
			case 't':
				threshold = atof(optarg);
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

		quit = true;

		cerr << "Usage for bigsi++ (v. " << BIGSI_VERSION << "):" << endl;
		cerr << "\t[-o <output file>] (default is stdout)" << endl;
		cerr << "\t[--o.csv (output CSV) | --o.json (output JSON)]" << endl;
		cerr << "\t[-t <search threshold>] (default is " 
			<< DEFAULT_SEARCH_THRESHOLD << ")" << endl;
		cerr << "\t-d <database search path> (can be repeated)" << endl;
		cerr << "\t[-i <input sequence file>] (can be repeated)" << endl;
		cerr << "\t[<DNA sequence>] (can be repeated)" << endl;
		return;
	}

	for(int i = optind;i < argc;i++){
		query_seq.push_back(argv[i]);
	}
	
	// Recursively add any database files in the provided search path(s)
	while( ff.next() ){

		if( find_file_extension(ff.name(), allowed_db_extentions) ) {
			subject_files.push_back( ff.name() );
		}
	}

	if( subject_files.empty() ){
		
		cerr << "Please provide at least one database file to search (-d)" << endl;
		quit = true;
		return;
	}
	
	if( query_files.empty() && query_seq.empty() ){
		
		cerr << "Please provide at least one query sequence or file" << endl;
		quit = true;
		return;
	}
	
	if( (threshold <= 0.0) || (threshold > 1.0) ){
		
		cerr << "Please provide: 0.0 < search threshold <= 1.0" << endl;
		quit = true;
		return;
	}
}

bool is_fasta(const string &m_input)
{

	if( match_extension(m_input, ".fna") || match_extension(m_input, ".fna.gz") ){
		return true;
	}

	if( match_extension(m_input, ".fa") || match_extension(m_input, ".fa.gz") ){
                return true;
        }

	if( match_extension(m_input, ".fasta") || match_extension(m_input, ".fasta.gz") ){
                return true;
        }

	return false;
}

size_t str_to_uint32_t(const string &m_str)
{
	size_t ret = 0;
	size_t power = 1;

	for(string::const_reverse_iterator i = m_str.rbegin();i != m_str.rend();++i){
	
		if( (*i < '0') || (*i > '9') ){
			throw __FILE__ ":str_to_uint32_t: Illegal character";
		}
	
		ret += power*(*i - '0');
		power *= 10;		
	}

	if(ret > UINT_MAX){
		throw __FILE__ ":str_to_uint32_t: Overflow!";
	}
	
	return uint32_t(ret);

}

void DownloadOptions::load(int argc, char* argv[])
{	
	// Don't quit unless we need to
	quit = false;
	
	list_only = false;
	false_positive_probability = DEFAULT_FALSE_POSITIVE_PROBABILITY;
	kmer_len = DEFAULT_KMER_LENGTH;
	min_kmer_count = DEFAULT_SRA_MIN_KMER_COUNT;
	max_num_download_attempts = DEFAULT_DOWNLOAD_ATTEMPT;
	hash_func = MURMUR_HASH;
	sleep_interval = 0;
	max_backlog = DEFAULT_MAX_BACKLOG;
	num_download_threads = DEFAULT_DOWNLOAD_THREADS;
	
	const char* options = "i:k:p:t:h?";
	int config_opt = 0;
	int long_index = 0;

	struct option long_opts[] = {
		{"download", true, &config_opt, 1},
		{"bloom", true, &config_opt, 2},
		{"list", false, &config_opt, 3},
		{"min-kmer-count", true, &config_opt, 4},
		{"log", true, &config_opt, 5},
		{"max-retry", true, &config_opt, 6},
		{"hash", true, &config_opt, 7},
		{"sleep", true, &config_opt, 9},
		{"max-backlog", true, &config_opt, 10},
		{0,0,0,0} // Terminate options list
	};
	
	int opt_code;
	opterr = 0;

	bool print_usage = (argc == 1);

	while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

		switch( opt_code ){
			case 0:
				if(config_opt == 1){ // --download
						
					download_dir = optarg;
					break;
				}
				
				if(config_opt == 2){ // --bloom
						
					bloom_dir = optarg;
					break;
				}
							
				if(config_opt == 3){ // --list
						
					list_only = true;
					break;
				}
				
				if(config_opt == 4){ // --min-kmer-count
						
					min_kmer_count = str_to_uint32_t(optarg);
					break;
				}

				if(config_opt == 5){ // --log
						
					log_file = optarg;
					break;
				}

				if(config_opt == 6){ // --max-retry
						
					max_num_download_attempts = str_to_uint32_t(optarg);
					break;
				}

				if(config_opt == 7){ // --hash
						
					hash_func = parse_hash_function_name(optarg);
					break;
				}
								
				if(config_opt == 9){ // --sleep
						
					sleep_interval = str_to_uint32_t(optarg);
					break;
				}
				
				if(config_opt == 10){ // --max-backlog
						
					max_backlog = str_to_uint32_t(optarg);
					break;
				}
				
				cerr << "Unknown flag!" << endl;
				break;
			case 'i':
				metadata_file = optarg;
				break;
			case 'k':
				kmer_len = str_to_uint32_t(optarg);
				break;
			case 'p':
				false_positive_probability = atof(optarg);
				break;
			case 't':
				num_download_threads = str_to_uint32_t(optarg);
				break;
			case '?':
			case 'h':
				print_usage = true;
				break;
			default:
				cerr << '\"' << (char)opt_code << "\" is not a valid option!" << endl;
				break;
		};
	}

	if(print_usage){

		quit = true;

		cerr << "Usage for SRA Download (v. " << DOWNLOAD_VERSION << "):" << endl;
		cerr << "\t-i <XML metadata file from ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata>" << endl;
		cerr << "\t--download <download directory for SRA data>" << endl;
		cerr << "\t--bloom <Bloom filter output directory>" << endl;
		cerr << "\t--log <log filename> (required to track progress for restarts)" << endl;
		cerr << "\t[-k <kmer length>] (default is " << DEFAULT_KMER_LENGTH << ")" << endl;
		cerr << "\t[-p <false positive probability (per k-mer, per-filter)>] (default is "
			<< DEFAULT_FALSE_POSITIVE_PROBABILITY << ")" << endl;
		cerr << "\t[--min-kmer-count <minimum allowed k-mer count>] (default is " 
			<< DEFAULT_SRA_MIN_KMER_COUNT << ")" << endl;
		cerr << "\t[--hash <hash function name>] (default is " 
			<< hash_name(hash_func) << ")" << endl;
		
		cerr << "\t\tAllowed hash functions: ";
		
		for(int i = 0;i < UNKNOWN_HASH;++i){
			
			if(i != 0){
				cerr << ", ";
			}
			
			cerr << hash_name( HashFunction(i) );
		}
		
		cerr << endl;
		
		cerr << "\t[--list (list, but do not download, SRA data)]" << endl;
		cerr << "\t[--sleep <sec> (time to sleep between downloads)]" << endl;
		cerr << "\t[--max-backlog <number of downloads> (Pause SRA downloading when exceeded)]" 
			" (default is " << DEFAULT_MAX_BACKLOG << "; 0 is no limit)" << endl;
		cerr << "\t[-t <number of download threads>] (default is " 
			<< DEFAULT_DOWNLOAD_THREADS<< ")" << endl;
		
		return;
	}

	if(kmer_len == 0){

		quit = true;
		cerr << "Please enter a kmer size (-k) > 0" << endl;
		return;
	}
	
	if(kmer_len > MAX_WORD_LEN){

		quit = true;
		cerr << "Please enter a kmer size (-k) < " << MAX_WORD_LEN << endl;
		return;
	}
	
	if(num_download_threads == 0){
		
		quit = true;
		cerr << "Please specify at least one thread to use for downloading SRA data(-t)" 
			<< endl;
		return;
	}
	
	if( (false_positive_probability <= 0.0f) || (false_positive_probability >= 1.0f) ){
	
		quit = true;
		cerr << "Please enter a per-kmer, per-filter false positive probability: 0.0 < p < 1.0" << endl;
		return;
	}

	if( log_file.empty() ){

		quit = true;
		cerr << "Please specify a log file (--log)" << endl;
		return;
	}
	
	if( metadata_file.empty() ){

		quit = true;
		cerr << "Please specify a metadata input file (-i)" << endl;
		return;
	}

	if( !list_only && download_dir.empty() ){
		
		quit = true;
		cerr << "Please provide a download directory (--download)" << endl;
		return;
	}
	
	if( !list_only && bloom_dir.empty() ){
		
		quit = true;
		cerr << "Please provide a Bloom filter directory (--bloom)" << endl;
		return;
	}
	
	if( hash_func == UNKNOWN_HASH ){
		
		quit = true;
		
		cerr << "Please specify a valid hash function (one of: ";
		
		for(int i = 0;i < UNKNOWN_HASH;++i){
			
			if(i != 0){
				cerr << ", ";
			}
			
			cerr << hash_name( HashFunction(i) );
		}
		
		cerr << ")";
		
		return;
	}
}

void BloomerOptions::load(int argc, char* argv[])
{	
	// Don't quit unless we need to
	quit = false;
	
	false_positive_probability = DEFAULT_FALSE_POSITIVE_PROBABILITY;
	kmer_len = DEFAULT_KMER_LENGTH;
	min_kmer_count = DEFAULT_SRA_MIN_KMER_COUNT;
	num_file_slice = DEFAULT_NUM_SRA_FILE_SLICE;
	hash_func = MURMUR_HASH;
	save_sra = false; // Should we remove the SRA data after filter construction?
	verbose = false;
	read_meta_data = true;
	
	const char* options = "i:o:k:p:vh?";
	int config_opt = 0;
	int long_index = 0;

	struct option long_opts[] = {
		{"min-kmer-count", true, &config_opt, 1},
		{"hash", true, &config_opt, 2},
		{"save-sra", false, &config_opt, 3},
		{"slice", true, &config_opt, 4},
		{"no-meta", false, &config_opt, 5},
		{0,0,0,0} // Terminate options list
	};
	
	int opt_code;
	opterr = 0;

	bool print_usage = (argc == 1);

	while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

		switch( opt_code ){
			case 0:
				
				if(config_opt == 1){ // --min-kmer-count
						
					min_kmer_count = str_to_uint32_t(optarg);
					break;
				}

				if(config_opt == 2){ // --hash
						
					hash_func = parse_hash_function_name(optarg);
					break;
				}
				
				if(config_opt == 3){ // --save-sra
						
					save_sra = true;
					break;
				}
				
				if(config_opt == 4){ // --slice
						
					num_file_slice = str_to_uint32_t(optarg);
					break;
				}
				
				if(config_opt == 5){ // --no-meta
						
					read_meta_data = false;
					break;
				}
				
				cerr << "Unknown flag!" << endl;
				break;
			case 'i':
				input_dir = optarg;
				break;
			case 'o':
				output_file = optarg;
				break;
			case 'k':
				kmer_len = str_to_uint32_t(optarg);
				break;
			case 'p':
				false_positive_probability = atof(optarg);
				break;
			case 'v':
				verbose = true;
				break;
			case '?':
			case 'h':
				print_usage = true;
				break;
			default:
				cerr << '\"' << (char)opt_code << "\" is not a valid option!" << endl;
				break;
		};
	}

	if(print_usage){

		quit = true;

		cerr << "Usage for Bloomer (v. " << BIGSI_VERSION << "):" << endl;
		cerr << "\t-i <input SRA directory>" << endl;
		cerr << "\t-o <Bloom filter output file>" << endl;
		cerr << "\t[-k <kmer length>] (default is " << DEFAULT_KMER_LENGTH << ")" << endl;
		cerr << "\t[-p <false positive probability (per k-mer, per-filter)>] (default is "
			<< DEFAULT_FALSE_POSITIVE_PROBABILITY << ")" << endl;
		cerr << "\t[--min-kmer-count <minimum allowed k-mer count>] (default is " 
			<< DEFAULT_SRA_MIN_KMER_COUNT << ")" << endl;
		cerr << "\t[--hash <hash function name>] (default is " 
			<< hash_name(hash_func) << ")" << endl;
		
		cerr << "\t\tAllowed hash functions: ";
		
		for(int i = 0;i < UNKNOWN_HASH;++i){
			
			if(i != 0){
				cerr << ", ";
			}
			
			cerr << hash_name( HashFunction(i) );
		}
		
		cerr << endl;
		
		//cerr << "\t[--save-sra (do not delete the input SRA directory)]" << endl;
		cerr << "\t[--no-meta (do not attempt to read the SRA meta-data file)]" << endl;
		cerr << "\t[--slice <num slice>] (number of SRA file slices to read in parallel; default is "
			<< DEFAULT_NUM_SRA_FILE_SLICE << ")" << endl;
		cerr << "\t[-v (turn on verbose output)]" << endl;
		
		return;
	}

	if( input_dir.empty() ){
		
		quit = true;
		cerr << "Please specify an input directory (-i)" << endl;
		return;
	}
	
	if( output_file.empty() ){
		
		quit = true;
		cerr << "Please specify an output Bloom filter file (-o)" << endl;
		return;
	}
	
	if(kmer_len == 0){

		quit = true;
		cerr << "Please enter a kmer size (-k) > 0" << endl;
		return;
	}
	
	if(kmer_len > MAX_WORD_LEN){

		quit = true;
		cerr << "Please enter a kmer size (-k) < " << MAX_WORD_LEN << endl;
		return;
	}
	
	if( (num_file_slice == 0) || (num_file_slice > MAX_SRA_FILE_SLICE) ){
		
		quit = true;
		cerr << "The number of SRA file slices (" << num_file_slice 
			<< ") is out of bounds 0 < num slice <= " 
			<< MAX_SRA_FILE_SLICE << endl;
		return;
	}
	
	if( (false_positive_probability <= 0.0f) || (false_positive_probability >= 1.0f) ){
	
		quit = true;
		cerr << "Please enter a per-kmer, per-filter false positive probability: 0.0 < p < 1.0" << endl;
		return;
	}

	if( hash_func == UNKNOWN_HASH ){
		
		quit = true;
		
		cerr << "Please specify a valid hash function (one of: ";
		
		for(int i = 0;i < UNKNOWN_HASH;++i){
			
			if(i != 0){
				cerr << ", ";
			}
			
			cerr << hash_name( HashFunction(i) );
		}
		
		cerr << ")";
		
		return;
	}
}
