// Bloom filter test rig
#include <iostream>
#include <getopt.h>
#include <math.h>

#include "maestro.h"
#include "string_conversion.h"

using namespace std;

// Keep MPI happy
int mpi_rank;
int mpi_numtasks;

int main(int argc, char *argv[])
{
	try{

		const char* options = "d:k:p:h?";
		int config_opt = 0;
		int long_index = 0;

		struct option long_opts[] = {
			{"min-kmer-count", true, &config_opt, 1},
			{"hash", true, &config_opt, 2},
			{"len.min", true, &config_opt, 6},
			{"len.max", true, &config_opt, 7},
			{0,0,0,0} // Terminate options list
		};

		int opt_code;
		opterr = 0;

		bool print_usage = (argc == 1);

		string output_dir = ".";
		deque<string> sra_accessions;
		MaestroOptions opt;

		opt.status_file = DEFAULT_STATUS_FILE;	
		opt.false_positive_probability = DEFAULT_FALSE_POSITIVE_PROBABILITY;
		opt.kmer_len = DEFAULT_KMER_LENGTH;
		opt.min_kmer_count = DEFAULT_SRA_MIN_KMER_COUNT;
		opt.min_log_2_filter_len = DEFAULT_MIN_LOG_2_FILTER_LEN;
		opt.max_log_2_filter_len = DEFAULT_MAX_LOG_2_FILTER_LEN;
		opt.hash_func = MURMUR_HASH_32;
		opt.num_download_attempt = DEFAULT_DOWNLOAD_ATTEMPT;
		opt.limit_num_download = 0; // <-- By default, *don't* limit
		opt.max_sra_file_size_GB = DEFAULT_MAX_SRA_FILE_SIZE;
		opt.s3_no_write = false; // By default, we write the database files to S3
		opt.save_bloom = false;
		opt.save_db = false;
		opt.save_sra = false;
		opt.stream_sra = false; // By default, first download using prefetch
		opt.verbose = false;
		
		while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

			switch( opt_code ){
				case 0:
					
					if(config_opt == 1){ // --min-kmer-count
						
						opt.min_kmer_count = str_to_uint32_t(optarg);
						break;
					}

					if(config_opt == 2){ // --hash
							
						opt.hash_func = parse_hash_function_name(optarg);
						break;
					}

					if(config_opt == 6){ // --len.min
						
						opt.min_log_2_filter_len = str_to_uint32_t(optarg);
						break;
					}
					
					if(config_opt == 7){ // --len.max
							
						opt.max_log_2_filter_len = str_to_uint32_t(optarg);
						break;
					}
					cerr << "Unknown flag!" << endl;
					break;
				case 'd':
					output_dir = optarg;
					break;
				case 'k':
					opt.kmer_len = str_to_uint32_t(optarg);
					break;
				case 'p':
					opt.false_positive_probability = atof(optarg);
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

		for(int i = optind;i < argc;i++){
			sra_accessions.push_back(argv[i]);
		}

		if(print_usage){

			cerr << "Usage for Bloom filter test rig:" << endl;
			cerr << "\t[-d <output directory for Bloom filter>] (default is '.')" << endl;
			cerr << "\t[-k <kmer length>] (default is " << DEFAULT_KMER_LENGTH << ")" << endl;
			cerr << "\t[-p <false positive probability (per k-mer, per-filter)>] (default is "
				<< DEFAULT_FALSE_POSITIVE_PROBABILITY << ")" << endl;
			cerr << "\t[--min-kmer-count <minimum allowed k-mer count>] (default is " 
				<< DEFAULT_SRA_MIN_KMER_COUNT << ")" << endl;
			cerr << "\t[--hash <hash function name>] (default is " 
				<< hash_name(opt.hash_func) << ")" << endl;
				
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
			cerr << "\t<SRA accession to read1> ..." << endl;
			return EXIT_SUCCESS;
		}

		if( sra_accessions.empty() ){

			cerr << "Please specify one or more SRA accessions to read" << endl;
			return EXIT_SUCCESS;
		}

		cerr << "kmer len = " << opt.kmer_len << endl;
		cerr << "min kmer count = " << opt.min_kmer_count << endl;
		cerr << "min log2 filter len = " << opt.min_log_2_filter_len << endl;
		cerr << "max log2 filter len = " << opt.max_log_2_filter_len << endl;
		cerr << "false positive probability = " << opt.false_positive_probability << endl;
		
		for(deque<string>::const_iterator acc = sra_accessions.begin();acc != sra_accessions.end();++acc){

			time_t profile = time(NULL);

			BloomParam param;
			FilterInfo info;
			BloomProgress progress;

			cerr << "Building Bloom filter for " << *acc << endl;

			const uint64_t num_bp = number_of_bases(*acc);

			cerr << "Number of bases (from metadata) = " << num_bp << endl;

			if(num_bp > 0){

				cerr << "log(bp) = " << log(num_bp)/log(2.0) << endl;

				// The desired false positive rate for the counting filter
				//const double fp_count = 1.0e-3;
				const double fp_count = 1.0e-2;

				// Required Counting Bloom filter length assuming two counting filters, each with
				// two hash functions (and the same number of bits).
				double counting_length = 1.0/( 1.0 - pow( 1.0 - pow(fp_count, 1.0/4.0) , 1.0/(2*num_bp) ) );

				cerr << "Desired counting length = " << counting_length << "; log(L) = " 
					<< log(counting_length)/log(2.0) << endl;
			}

			unsigned char ret = make_bloom_filter(str_to_accession(*acc), 
				info, param, progress, output_dir, opt);

			// Handle aligned colorspace reads as a special case. According to https://github.com/ncbi/ncbi-vdb/issues/31
			// aligned colorspace SRA records are "broken" and will fail when trying to read primary alignments 
			// followed by unaligned reads. The signature of this failure is that all primary alignments will be
			// read successfully and no unaligned reads will be successfully read.
			if( (ret == STATUS_BLOOM_FAIL) && (progress.num_primary_align > 0) && 
				(progress.curr_primary_align == progress.num_primary_align) &&
				(progress.num_unaligned_read > 0) && (progress.curr_unaligned_read == 0) ){
				
				cerr << "\tColor space edge case; forcing unaligned read loading" << endl;

				ret = make_bloom_filter(str_to_accession(*acc), info, param, 
					progress, output_dir, opt, true /*force unaligned*/);
			}

			profile = time(NULL) - profile;

			cerr << "\tmake_bloom_filter(" << *acc << ") returned \"";

			switch(ret){
				case STATUS_BLOOM_SUCCESS:
					cerr << "success\" in " << profile << " sec" << endl;
					break;
				case STATUS_BLOOM_FAIL:
					cerr << "failure\" in " << profile << " sec" << endl;
					break;
				case STATUS_BLOOM_INVALID:
					cerr << "invalid\" in " << profile << " sec" << endl;
					break;
				default:
					cerr << "unknown\" in " << profile << " sec" << endl;
			};

			cerr << "\tnum_bp = " << progress.num_bp << endl;
			cerr << "\tnum_kmer = " << progress.num_kmer << endl;
			cerr << "\tlog2 filter length = " << param.log_2_filter_len << endl;
			cerr << "\tnum hash functions = " << param.num_hash << endl;

			cerr << "\tvalid_read_collection = " << (progress.valid_read_collection ? "true" : "false") << endl;

			if(progress.num_primary_align > 0){

				cerr << "\tPrimary alignment: " << progress.curr_primary_align << " out of " 
					<< progress.num_primary_align << " (" 
					<< (100.0*progress.curr_primary_align)/progress.num_primary_align << "%)" << endl;
			}

			if(progress.num_unaligned_read > 0){

				cerr << "\tUnaligned reads: " << progress.curr_unaligned_read << " out of " 
					<< progress.num_unaligned_read << " (" 
					<< (100.0*progress.curr_unaligned_read)/progress.num_unaligned_read << "%)" << endl;

				cerr << "\tCurrent fragment: " << progress.curr_fragment << endl;
			}

			if(progress.num_read > 0){

				cerr << "\tReads: " << progress.curr_read << " out of " 
					<< progress.num_read << " (" << (100.0*progress.curr_read)/progress.num_read << "%)" << endl;
				cerr << "\tCurrent fragment: " << progress.curr_fragment << endl;
			}

			if( !progress.error.empty() ){
				cerr << "\tError: " << progress.error << endl;
			}
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