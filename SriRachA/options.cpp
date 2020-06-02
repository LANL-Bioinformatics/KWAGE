#include <iostream>

#include <getopt.h>
#include <math.h>

#include "options.h"
#include "sriracha.h"

using namespace std;

void SrirachaOptions::load(int argc, char* argv[])
{
	// Don't quit unless we need to
	quit = false;
	
	sra_accession.clear();
	input_sequence_files.clear();
	output_filename.clear();
	sra_accession_filename.clear();
	search_strategy = SEARCH_BY_KMER;
	kmer_len = DEFAULT_KMER_LENGTH;
	kmer_match_threshold = DEFAULT_KMER_MATCH_THRESHOLD;
	min_valid_kmer = DEFAULT_MIN_VALID_KMER;
	min_read_complexity = DEFAULT_MIN_READ_COMPLEXITY;
	min_read_length = DEFAULT_MIN_READ_LENGTH;
	max_num_match = DEFAULT_MAX_MATCH;
	max_retry = 0;

	// By default, read the entire SRA record
	slice_index = 0;
	num_slice = 1;

	verbose = SILENT;

	// Command line arguments
	// -i <input sequence files> (can be repeated)
	// [-o <output filename>] (default is stdout)
	// [-a <list of SRA accessions in a text file>]
	// [--read.len.min <minimum read length>] (default is DEFAULT_MIN_READ_LENGTH)
	// [-v] (increase the verbosity; default is silent)
	// [--max-results <maximum number of results to show per accession/target>] (default is DEFAULT_MAX_MATCH)
	// [--retry <maximum number of download atttemps>] (default is 0)
	// [--slice <slice number [0, N)]>] (no MPI allowed!)
	// [--of <number of slices, N>] (no MPI allowed!)
	// Search strategy
	// [--search-by-align]
	// [--search-by-kmer] (default)
	// 		[-k <k-mer length>] (default is DEFAULT_KMER_LENGTH)
	//		[-t <match threshold>] (default is DEFAULT_KMER_MATCH_THRESHOLD)
	//		[-n <min number valid kmer>] (default is DEFAULT_MIN_VALID_KMER)]
	//		[--read.complexity.min <min read complexity>] (default is DEFAULT_MIN_READ_COMPLEXITY)
	// [--search-by-bloom]
	// <SRA accession1> ...

	const char* options = "k:t:n:o:i:a:vh?";
	int config_opt = 0;
	int long_index = 0;

	struct option long_opts[] = {
		{"search-by-align", false, &config_opt, 1},
		{"search-by-kmer", false, &config_opt, 2},
		{"search-by-bloom", false, &config_opt, 3},
		{"read.len.min", true, &config_opt, 4},
		{"read.complexity.min", true, &config_opt, 5},
		{"max-results", true, &config_opt, 6},
		{"vv", false, &config_opt, 7},
		{"vvv", false, &config_opt, 8},
		{"vvvv", false, &config_opt, 9},
		{"retry", true, &config_opt, 10},
		{"slice", true, &config_opt, 11},
		{"of", true, &config_opt, 12},
		{0,0,0,0} // Terminate options list
	};

	int opt_code;
	opterr = 0;

	bool print_usage = (argc == 1);

	while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

		switch( opt_code ){
			case 0:
				
				if(config_opt == 1){ // search-by-align

					search_strategy = SEARCH_BY_ALIGN;
					break;
				}

				if(config_opt == 2){ // search-by-kmer

					search_strategy = SEARCH_BY_KMER;
					break;
				}

				if(config_opt == 3){ // search-by-bloom

					search_strategy = SEARCH_BY_BLOOM;
					break;
				}

				if(config_opt == 4){ // read.len.min

					min_read_length = abs( atoi(optarg) );
					break;
				}

				if(config_opt == 5){ // read.complexity.min

					min_read_complexity = atof(optarg);
					break;
				}

				if(config_opt == 6){ // max-results

					// Read as float and cast to unsigned int
					max_num_match = fabs( atof(optarg) );
					break;
				}

				if(config_opt == 7){ // vv

					verbose += 2;
					break;
				}

				if(config_opt == 8){ // vvv

					verbose += 3;
					break;
				}

				if(config_opt == 9){ // vvvv

					verbose += 4;
					break;
				}

				if(config_opt == 10){ // retry

					max_retry = abs( atoi(optarg) );
					break;
				}

				if(config_opt == 11){ // slice

					slice_index = abs( atoi(optarg) );
					break;
				}

				if(config_opt == 12){ // of

					num_slice = abs( atoi(optarg) );
					break;
				}

				cerr << "Unknown flag!" << endl;
				break;
			case 'a':
				// Allow the user to specify a test file containing a 
				// list of SRA accessions to search
				sra_accession_filename = optarg;
				break;
			case 'n':
				min_valid_kmer = abs( atoi(optarg) );
				break;
			case 't':
				kmer_match_threshold = atof(optarg);
				break;
			case 'k':
				kmer_len = abs( atoi(optarg) );
				break;
			case 'i':
				input_sequence_files.push_back(optarg);
				break;
			case 'o':
				output_filename = optarg;
				break;
			case '?':
			case 'h':
				print_usage = true;
				break;
			case 'v':
				++verbose;
				break;
			default:
				cerr << '\"' << (char)opt_code << "\" is not a valid option!" << endl;
				break;
		};
	}

	for(int i = optind;i < argc;i++){
		sra_accession.push_back(argv[i]);
	}

	if(print_usage){

		quit = true;

		cerr << "Usage for SriRachA (v. " << SRIRACHA_VERSION << "):" << endl;
		cerr << "\t-i <input sequence files> (can be repeated)" << endl;
		cerr << "\t[-o <output filename>] (default is stdout)" << endl;
		cerr << "\t[--read.len.min <minimum read length>] (default is " << DEFAULT_MIN_READ_LENGTH << ")" << endl;
		cerr << "\t[--max-results <maximum number of results to show per accession/query>] (default is " 
			<< DEFAULT_MAX_MATCH << ")" << endl;
		cerr << "\t[-a <list of SRA accessions in a text file>]" << endl;
		cerr << "\t[-v (increase the verbosity: silent, tacitern, normal, chatty. Default is silent)]" << endl;
		cerr << "\t[--retry <maximum number of download atttemps>] (default is 0)" << endl;
		cerr << "\t[--slice <slice number [0, N)]>] (not compatible with MPI)" << endl;
		cerr << "\t[--of <number of slices, N>] (not compatible with MPI)" << endl;
		cerr << "\tSearch strategies" << endl;
		cerr << "\t\t[--search-by-align]" << endl;
		cerr << "\t\t\t(Not implemented yet!)" << endl;
		cerr << "\t\t[--search-by-kmer] (default)" << endl;
		cerr << "\t\t\t[-k <k-mer length>] (default is " << DEFAULT_KMER_LENGTH << ")" << endl;
		cerr << "\t\t\t[-t <match threshold>] (default is " << DEFAULT_KMER_MATCH_THRESHOLD << ")" << endl;
		cerr << "\t\t\t[-n <min number valid kmer>] (default is " << DEFAULT_MIN_VALID_KMER << ")" << endl;
		cerr << "\t\t\t[--read.complexity.min <min read complexity>] (default is " << DEFAULT_MIN_READ_COMPLEXITY << ")" << endl;
		cerr << "\t\t[--search-by-bloom]" << endl;
		cerr << "\t\t\t(Not implemented yet!)" << endl;
		cerr << "\t<SRA accession or file or dir> ..." << endl;
		
		return;
	}

	if(min_valid_kmer == 0){

		quit = true;

		cerr << "Please specify: 0 < minimum number of valid kmers" << endl;
		
		return;
	}

	if(max_num_match == 0){

		quit = true;

		cerr << "Please specify: 0 < max number of matches to report" << endl;
		
		return;
	}

	if( (kmer_len < MIN_KMER_LEN) || (kmer_len > MAX_KMER_LEN) ){

		quit = true;

		cerr << "Please specify: " << MIN_KMER_LEN 
			<< " <= kmer length <= " << MAX_KMER_LEN << endl;
		
		return;
	}

	if( (kmer_match_threshold <= 0.0f) || (kmer_match_threshold > 1.0f) ){

		quit = true;

		cerr << "Please specify: 0.0 < kmer match threshold <= 1.0" << endl;
		
		return;
	}

	if( (min_read_complexity <= 0.0f) || (min_read_complexity > 1.0f) ){

		quit = true;

		cerr << "Please specify: 0.0 < min read complexity <= 1.0" << endl;
		
		return;
	}

	if( input_sequence_files.empty() ){

		quit = true;

		cerr << "Please specify one or more input sequence files (-i)" << endl;
		
		return;
	}

	// If the user does not specify an SRA accession, we assume that they would like
	// to pipe from stdin
	//if( sra_accession.empty() ){
	//	quit = true;
	//	cerr << "Please specify one or more SRA run accessions" << endl;	
	//	return;
	//}

	if(search_strategy == UNDEFINED_SEARCH){

		quit = true;

		cerr << "Please specify a valid search strategy" << endl;
		
		return;
	}

	if(slice_index >= num_slice){

		quit = true;

		cerr << "Please specify slice index (--slice) less than the number of slices (--of)" << endl;
		
		return;
	}
}