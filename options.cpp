#include <iostream> 
#include <unordered_set>
#include <unordered_map>
#include <deque>
#include <algorithm>

#include <getopt.h>
#include <string.h>
#include <math.h>

#include "options.h"
#include "word.h"
#include "file_util.h"
#include "kwage.h"
#include "hash.h"
#include "string_conversion.h"
#include "maestro.h"

using namespace std;

// Enumerate the allowed file extensions for the KWAGE query sequence(s)
const char* allowed_sequence_extentions [] = {
	".fna", ".fna.gz",
	".fasta", ".fasta.gz",
	".fa", ".fa.gz",
	".fastq", ".fastq.gz",
	NULL
};

const char* allowed_db_extentions [] = {
	".db",
	NULL
};

const char* allowed_filter_extention = ".bloom";

bool is_fasta(const string &m_input);

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

		cerr << "Usage for KWAGE (v. " << KWAGE_VERSION << "):" << endl;
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
	
	for(deque<string>::const_iterator i = query_files.begin();i != query_files.end();++i){

		// Make sure that the specified filename matches one of the allow query sequence file
		// extensions
		bool valid = false;

		for(const char** ext = allowed_sequence_extentions;*ext != NULL;++ext){

			const size_t ext_len = strlen(*ext);

			if( ext_len > i->size() ){
				continue;
			}

			if( i->find(*ext) == (i->size() - ext_len) ){

				valid = true;
				break;
			}
		}

		if(!valid){
			
			cerr << "The query sequence file name, " << *i 
				<< ", does not have an allowed file extension" << endl;

			quit = true;
			return;
		}
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

void InventoryOptions::load(int argc, char* argv[])
{	
	// Don't quit unless we need to
	quit = false;
	
	metadata_file.clear();
	output_file.clear();
	required_strategy.clear();
	required_source.clear();
	include_accessions.clear();

	list_only = false;

	// By default, there is no date restriction. Make the begin
	// and end date far in the past and far in the future respectively
	begin_date = Date("0001-01-01"); // Jan 1, AD 1
	end_date = Date("9999-01-01"); // Jan 1, AD 9999

	const char* options = "i:o:h?";
	int config_opt = 0;
	int long_index = 0;

	struct option long_opts[] = {
		{"list", false, &config_opt, 3},
		{"date.from", true, &config_opt, 13},
		{"date.to", true, &config_opt, 14},
		{"strategy", true, &config_opt, 15},
		{"source", true, &config_opt, 16},
		{"include", true, &config_opt, 17},
		{0,0,0,0} // Terminate options list
	};
	
	int opt_code;
	opterr = 0;

	bool print_usage = (argc == 1);

	string date_from;
	string date_to;
	string include_filename; // An optional list of SRA accessions to include

	while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

		switch( opt_code ){
			case 0:
				if(config_opt == 3){ // --list
						
					list_only = true;
					break;
				}
				
				if(config_opt == 13){ // --date.from
						
					date_from = optarg;
					break;
				}

				if(config_opt == 14){ // --date.to
						
					date_to = optarg;
					break;
				}

				if(config_opt == 15){ // --strategy
						
					required_strategy.insert(optarg);
					break;
				}

				if(config_opt == 16){ // --source
						
					required_source.insert(optarg);
					break;
				}

				if(config_opt == 17){ // --include
						
					include_filename = optarg;
					break;
				}

				cerr << "Unknown flag!" << endl;
				break;
			case 'i':
				metadata_file = optarg;
				break;
			case 'o':
				output_file = optarg;
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

		cerr << "Usage for SRA inventory (v. " << INVENTORY_VERSION << "):" << endl;
		cerr << "\t-i <XML metadata file from ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata>" << endl;
		cerr << "\t[-o <binary output file>]" << endl;			
		cerr << "\t[--list (list, but do not write binary SRA inventory)]" << endl;
		cerr << "\t[--date.from <YYYY-MM-DD>] (only download SRA records received after this date)" << endl;
		cerr << "\t[--date.to <YYYY-MM-DD>] (only download SRA records received before this date)" << endl;
		cerr << "\t[--strategy <strategy key word>] (only download SRA records that match one of the specified experimental strategies)" << endl;
		cerr << "\t\tExamples include: RNA-Seq, WGS, AMPLICON, Bisulfite-Seq, ... (case sensitive!)" << endl;
		cerr << "\t[--source <source key word>] (only download SRA records that match one of the specified exterimental sources)" << endl;
		cerr << "\t\tExamples include: TRANSCRIPTOMIC, GENOMIC, METAGENOMIC, METATRANSCRIPTOMIC, ... (case sensitive!)" << endl;
		cerr << "\t[--include <list of SRA run accessions>] (only download SRA records that match one of the specified SRA runs)" << endl;
		return;
	}
	
	if( metadata_file.empty() ){

		quit = true;
		cerr << "Please specify a metadata input file (-i)" << endl;
		return;
	}

	if( !list_only && output_file.empty() ){
		
		quit = true;
		cerr << "Please provide an output file (-o)" << endl;
		return;
	}
	
	if( !date_from.empty() ){
		try{
			begin_date = Date(date_from);
		}
		catch(...){
			
			quit = true;

			cerr << "Please specify a valid data (YYYY-MM-DD) for --date.from" << endl;

			return;
		}
	}

	if( !date_to.empty() ){
		try{
			end_date = Date(date_to);
		}
		catch(...){
			
			quit = true;

			cerr << "Please specify a valid data (YYYY-MM-DD) for --date.to" << endl;

			return;
		}
	}

	if( !include_filename.empty() ){

		ifstream fin(include_filename);

		if(!fin){

			quit = true;

			cerr << "Unable to open the --include file to read a list of SRA run accessions" << endl;

			return;
		}

		string line;

		// Read a list of SRA runs accessions, one accession per line
		while( getline(fin, line) ){

			if( line.empty() ){
				continue;
			}

			include_accessions.push_back( str_to_accession(line) );
		}

		// Sort the list of accessions to include and make the list unqiue
		sort( include_accessions.begin(), include_accessions.end() );

		include_accessions.erase( unique( include_accessions.begin(), include_accessions.end() ), 
			include_accessions.end() );
	}
}

void MaestroOptions::load(int argc, char* argv[])
{	
	// Don't quit unless we need to
	quit = false;

	string scratch_dir;

	status_file = DEFAULT_STATUS_FILE;	
	false_positive_probability = DEFAULT_FALSE_POSITIVE_PROBABILITY;
	download_delay = DEFAULT_DOWNLOAD_DELAY;
	kmer_len = DEFAULT_KMER_LENGTH;
	min_kmer_count = DEFAULT_SRA_MIN_KMER_COUNT;
	min_log_2_filter_len = DEFAULT_MIN_LOG_2_FILTER_LEN;
	max_log_2_filter_len = DEFAULT_MAX_LOG_2_FILTER_LEN;
	hash_func = MURMUR_HASH_32;
	num_download_attempt = DEFAULT_DOWNLOAD_ATTEMPT;
	limit_num_download = 0; // <-- By default, *don't* limit
	max_sra_file_size_GB = DEFAULT_MAX_SRA_FILE_SIZE;
	retry_bloom = false;
	s3_no_write = false; // By default, we write the database files to S3
	save_bloom = false;
	save_db = false;
	save_sra = false;
	stream_sra = false; // By default, first download using prefetch
	verbose = false;

	const char* options = "k:p:vh?";
	int config_opt = 0;
	int long_index = 0;

	struct option long_opts[] = {
		{"min-kmer-count", true, &config_opt, 1},
		{"hash", true, &config_opt, 2},
		{"scratch", true, &config_opt, 3},
		{"s3", true, &config_opt, 4},
		{"meta", true, &config_opt, 5},
		{"len.min", true, &config_opt, 6},
		{"len.max", true, &config_opt, 7},
		{"status", true, &config_opt, 8},
		{"retry", true, &config_opt, 9},
		{"halt-after", true, &config_opt, 11},
		{"save.bloom", false, &config_opt, 12},
		{"save.db", false, &config_opt, 13},
		{"save.sra", false, &config_opt, 14},
		{"s3.no-write", false, &config_opt, 15},
		{"max-sra-download", true, &config_opt, 16},
		{"stream", false, &config_opt, 17},
		{"retry.bloom", false, &config_opt, 18},
		{"delay", true, &config_opt, 19},
		{"scratch.bloom", true, &config_opt, 20},
		{"scratch.database", true, &config_opt, 21},
		{"skip", true, &config_opt, 22},
		{0,0,0,0} // Terminate options list
	};
	
	int opt_code;
	opterr = 0;

	bool print_usage = (argc == 1);

	// Read the SRA accessions to skip as strings so we can check for
	// errors before converting to SraAccessions.
	deque<string> local_skip;

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

				if(config_opt == 3){ // --scratch
						
					scratch_dir = optarg;
					break;
				}

				if(config_opt == 4){ // --s3
						
					s3_bucket = optarg;
					break;
				}

				if(config_opt == 5){ // --meta
						
					metadata_file = optarg;
					break;
				}

				if(config_opt == 6){ // --len.min
						
					min_log_2_filter_len = str_to_uint32_t(optarg);
					break;
				}
				
				if(config_opt == 7){ // --len.max
						
					max_log_2_filter_len = str_to_uint32_t(optarg);
					break;
				}

				if(config_opt == 8){ // --status
						
					status_file = optarg;
					break;
				}

				if(config_opt == 9){ // --retry
						
					num_download_attempt = str_to_uint32_t(optarg);
					break;
				}

				if(config_opt == 11){ // --halt-after
						
					limit_num_download = str_to_uint32_t(optarg);
					break;
				}

				if(config_opt == 12){ // --save.bloom
						
					save_bloom = true;
					break;
				}

				if(config_opt == 13){ // --save.db
						
					save_db = true;
					break;
				}

				if(config_opt == 14){ // --save.sra
						
					save_sra = true;
					break;
				}

				if(config_opt == 15){ // --s3.no-write
						
					s3_no_write = true;
					break;
				}

				if(config_opt == 16){ // --max-sra-download
						
					max_sra_file_size_GB = str_to_uint32_t(optarg);
					break;
				}

				if(config_opt == 17){ // --stream
						
					stream_sra = true;
					break;
				}

				if(config_opt == 18){ // --retry.bloom
						
					retry_bloom = true;
					break;
				}

				if(config_opt == 19){ // --delay
						
					download_delay = atof(optarg);
					break;
				}

				if(config_opt == 20){ // --scratch.bloom
						
					scratch_bloom_dir = optarg;
					break;
				}

				if(config_opt == 22){ // --skip
						
					local_skip.push_back(optarg);
					break;
				}

				cerr << "Unknown flag!" << endl;
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

		cerr << "Usage for Maestro (v. " << MAESTRO_VERSION << "):" << endl;
		cerr << "** 32-bit hash functions limit the maximum Bloom filter length to 2^32 bits **" << endl;
		cerr << "\t--scratch <parent directory for staging Bloom filter and database files>" << endl;
		cerr << "\t[--scratch.bloom <scratch directory for staging Bloom filter>]" << endl;
		cerr << "\t[--scratch.database <scratch directory for staging database files>]" << endl;
		cerr << "\t--meta <binary metadata input file>" << endl;
		cerr << "\t--s3 <S3 path for storing database files>" << endl;
		cerr << "\t[--s3.no-write (do *not* write database files to s3)]" << endl;
		cerr << "\t[--stream (stream SRA data -- do not use prefetch to download!)]" << endl;
		cerr << "\t[--max-sra-download <max allowed SRA file size in GB>] (default is " << DEFAULT_MAX_SRA_FILE_SIZE
			<< " GB)" << endl;
		cerr << "\t[--status <binary SRA status file for restart>] (default is " << DEFAULT_STATUS_FILE << ")" << endl;
		cerr << "\t[--retry <number of download attempts>] (default is " << DEFAULT_DOWNLOAD_ATTEMPT << ")" << endl;
		cerr << "\t[--retry.bloom (retry all failed Bloom filters)]" << endl;
		cerr << "\t[--delay <minimum number of seconds between download/streaming requests>]" << endl;
		cerr << "\t[--halt-after <halt after this many SRA downloads> (default is not to stop)]" << endl;
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
		
		cerr << "\t[--len.min <log2 Bloom filter len>] (default is " 
			<< DEFAULT_MIN_LOG_2_FILTER_LEN << ")" << endl;
		cerr << "\t[--len.max <log2 Bloom filter len>] (default is " 
			<< DEFAULT_MAX_LOG_2_FILTER_LEN << ")" << endl;
			
		cerr << "\t[-v (turn on verbose output)]" << endl;
		cerr << "\t[--save.bloom (don't remove Bloom filters after database construction)]" << endl;
		cerr << "\t[--save.db (don't remove database file after S3 upload)]" << endl;
		cerr << "\t[--save.sra (don't remove SRA files after Bloom filter construction)]" << endl;
		cerr << "\t[--skip <SRA run accession> (skip over the specified accession; may be repeated)]" << endl;
		return;
	}
	
	if( metadata_file.empty() ){
		
		quit = true;
		cerr << "Please specify a binary metadata file of SRA accession information (--meta)" << endl;
		return;
	}

	if( scratch_dir.empty() ){
		
		if( scratch_bloom_dir.empty() && scratch_database_dir.empty() ){

			quit = true;
			cerr << "Please specify a scratch directory for staging Bloom filter and database files (--scratch)" << endl;
			return;
		}

		if( scratch_bloom_dir.empty() ){

			quit = true;
			cerr << "Please specify a scratch directory for staging Bloom filter files (--scratch.bloom)" << endl;
			return;
		}

		if( scratch_database_dir.empty() ){

			quit = true;
			cerr << "Please specify a scratch directory for staging database files (--scratch.database)" << endl;
			return;
		}
	}
	else{

		if( scratch_bloom_dir.empty() ){
			scratch_bloom_dir = scratch_dir + PATH_SEPARATOR + "bloom";
		}

		if( scratch_database_dir.empty() ){
			scratch_database_dir = scratch_dir + PATH_SEPARATOR + "database";
		}
	}

	if( s3_bucket.empty() ){
		
		quit = true;
		cerr << "Please specify the S3 bucket for storing database files (--s3)" << endl;
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
	
	if(min_kmer_count > MAX_SRA_MIN_KMER_COUNT){

		quit = true;
		cerr << "Please enter a min kmer count <= " << MAX_SRA_MIN_KMER_COUNT << " (--min-kmer-count)" << endl;
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
	
	if(min_log_2_filter_len > max_log_2_filter_len){
		
		quit = true;
		
		cerr << "Please specify a min log2 Bloom filter length less than the maximum filter length"
			<< endl;
		
		return;
	}
	
	if( (hash_func == MURMUR_HASH_32) && (max_log_2_filter_len > 32) ){
		
		quit = true;
		
		cerr << "When using a 32-bit hash function, please specify a maximum log2 Bloom filter length <= 32"
			<< endl;
		
		return;
	}

	// We currently only have status tags to allow upto 5 download attempts
	if(num_download_attempt > MAX_NUM_DOWNLOAD_FAIL){

		quit = true;
		
		cerr << "Please specify a number of download attempts <= " << MAX_NUM_DOWNLOAD_FAIL << endl;
		
		return;
	}

	if(max_sra_file_size_GB == 0){

		quit = true;
		
		cerr << "Please specify a maximum SRA file download size > 0 GB (--max-sra-download)" << endl;
		
		return;
	}

	for(deque<string>::const_iterator i = local_skip.begin();i != local_skip.end();++i){

		try{
			skip_sra.push_back( str_to_accession(*i) );
		}
		catch(...){

			cerr << "Unable to convert " << *i << " into a valid SRA run accession to skip" << endl;
			quit = true;
			return;
		}
	}

	sort( skip_sra.begin(), skip_sra.end() );

	// Remove duplicate accessions to skip
	skip_sra.erase( unique( skip_sra.begin(), skip_sra.end() ), skip_sra.end() );
}

ostream& operator<<(ostream &m_out, const MaestroOptions &m_opt)
{
	m_out << "Maestro v. " << MAESTRO_VERSION << endl;
	m_out << "Metadata input file = " << m_opt.metadata_file << endl;
	m_out << "Bloom scratch directory = " << m_opt.scratch_bloom_dir << endl;
	m_out << "Database scratch directory = " << m_opt.scratch_database_dir << endl;
	m_out << "Status file = " << m_opt.status_file << endl;
	m_out << "S3 bucket = " << m_opt.s3_bucket << endl;
	m_out << "Max SRA file size = " << m_opt.max_sra_file_size_GB << " GB" << endl;
	m_out << "False positive probability = " << m_opt.false_positive_probability << endl;
	m_out << "Minimum kmer count = " << m_opt.min_kmer_count << endl;
	m_out << "Kmer length = " << m_opt.kmer_len << endl;
	m_out << "Minimum Log_2 filter length = " << m_opt.min_log_2_filter_len << endl;
	m_out << "Maximu Log_2 filter length = " << m_opt.max_log_2_filter_len << endl;
	m_out << "Hash Function = " << hash_name(m_opt.hash_func) << endl;
	m_out << "Number of download attempts = " << m_opt.num_download_attempt << endl;

	if(m_opt.limit_num_download == 0){
		m_out << "Unlimited number of SRA downloads" << endl;
	}
	else{
		m_out << "Number of SRA downloads is limited to " << m_opt.limit_num_download << endl;
	}
	
	m_out << "Retrying previously failed Bloom filter creation = " << (m_opt.retry_bloom ? "true" : "false") << endl;
	m_out << "Save scratch copy of Bloom filter files = " << (m_opt.save_bloom ? "true" : "false") << endl;
	m_out << "Save scratch copy of database files = " << (m_opt.save_db ? "true" : "false") << endl;
	m_out << "Save scratch copy of SRA files = " << (m_opt.save_sra ? "true" : "false") << endl;
	m_out << "Write database files to S3 storage = " << (m_opt.s3_no_write ? "false": "true") << endl;
	m_out << "Streaming SRA data = " << (m_opt.stream_sra ? "true": "false") << endl;
	m_out << "Download/streaming delay = " << m_opt.download_delay << " sec" << endl;
	m_out << "Verbose output = " << (m_opt.verbose ? "true" : "false") << endl;

	return m_out;
}
