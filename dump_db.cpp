// Dump the contents of BIGSI++ database file in a human readable
// format.
// J. D. Gans
// Bioscience Division, B-11
// Los Alamos National Laboratory
// Mon Mar  2 12:02:07 2020

#include <iostream>
#include <iomanip>
#include <fstream>

#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include "bigsi++.h"
#include "bloom.h"
#include "string_conversion.h"
#include "keys.h"

using namespace std;

int main(int argc, char *argv[])
{
	try{
		
		// Command line arguments
		// -i <input database file> (can be repeated)
		// [-o <output file name>] (default is stdout)
		// [--bits <maximum number of bit-slices to display>]
		// [--bits.all] (display all bit-slices)
		// [--bits.none] (don't display any bit-slices) <-- default

		deque<string> input_filename;
		string output_filename;
		size_t num_bit_slice = 0;

		const char* options = "o:i:h?";
		int config_opt = 0;
		int long_index = 0;

		struct option long_opts[] = {
			{"bits", true, &config_opt, 1},
			{"bits.all", false, &config_opt, 2},
			{"bits.none", false, &config_opt, 3},
			{0,0,0,0} // Terminate options list
		};

		int opt_code;
		opterr = 0;

		bool print_usage = (argc == 1);
	
		while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

			switch( opt_code ){
				case 0:
					if(config_opt == 1){ // bits
							
						num_bit_slice = str_to_size_t(optarg);
						break;
					}
					
					if(config_opt == 2){ // bits.all
							
						num_bit_slice = 0xFFFFFFFFFFFFFFFF;
						break;
					}
					
					if(config_opt == 3){ // bits.none
							
						num_bit_slice = 0;
						break;
					}

					cerr << "Unknown flag!" << endl;
					break;
				case 'o':
					output_filename = optarg;
					break;
				case 'i':
					input_filename.push_back(optarg);
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
			cerr << "\t-i <input database file> (can be repeated)" << endl;
			cerr << "\t[-o <output file name>] (default is stdout)" << endl;
			cerr << "\t[--bits <maximum number of bit-slices to display>]" << endl;
			cerr << "\t[--bits.all] (display all bit-slices)" << endl;
			cerr << "\t[--bits.none] (don't display any bit-slices) <-- default" << endl;
			return EXIT_SUCCESS;
		}

		if( input_filename.empty() ){

			cerr << "Please specify one or more filenames" << endl;
			return EXIT_SUCCESS;
		}

		ofstream fout;

		if( !output_filename.empty() ){

			fout.open( output_filename.c_str() );

			if(!fout){

				cerr << "Unable to open output filename: " << output_filename << endl;
				return EXIT_FAILURE;
			}
		}

		ostream &out = fout.is_open() ? fout : cout;

		for(deque<string>::const_iterator f = input_filename.begin();f != input_filename.end();++f){
			
			ifstream fin(f->c_str(), ios::binary);
			
			if(!fin){
			
				cerr << "Unable to open " << *f << " for reading" << endl;
				return EXIT_FAILURE;
			}
			
			DBFileHeader header;
			
			binary_read(fin, header);
			
			if(!fin){
			
				cerr << "Unable to read database header" << endl;
				return EXIT_FAILURE;
			}
			
			out << "Header information for " << *f << endl;
			out << "\tmagic = " << header.magic << endl;
			out << "\tversion = " << header.version << endl;
			out << "\tcrc32 = " << std::hex << header.crc32 << std::dec << endl;
			out << "\tkmer_len = " << header.kmer_len << endl;
			out << "\tnum_hash = " << header.num_hash << endl;
			
			const size_t filter_len = header.filter_len();
			
			out << "\tfilter_len = " << filter_len << endl;
			
			out << "\tlog_2_filter_len = " << header.log_2_filter_len << endl;
			out << "\tnum_filter = " << header.num_filter << endl;
			
			switch(header.hash_func){
				case MURMUR_HASH_32:
					out << "\thash_func = Murmur32" << endl;
					break;
				case UNKNOWN_HASH:
					out << "\thash_func = Unknown" << endl;
					break;
				default:
					out << "\thash_func = Invalid" << endl;
					break;
			};
			
			switch(header.compression){
				case NO_COMPRESSION:
					out << "\tcompression = None" << endl;
					break;
				case RLE_COMPRESSION:
					out << "\tcompression = RLE" << endl;
					break;
				default:
					out << "\tcompression = Invalid" << endl;
					break;
			};
			
			if(header.compression != NO_COMPRESSION){

				cerr << "Compressed database files are not currently supported!" << endl;
				return EXIT_SUCCESS;
			}

			const unsigned int bytes_per_bitslice = header.num_filter/BitVector::BITS_PER_BLOCK +
				( (header.num_filter%BitVector::BITS_PER_BLOCK) > 0 ? 1 : 0);
			
			cout << "There are " << bytes_per_bitslice << " bytes per slice" << endl;
			
			cout << "Info start @ " << header.info_start << endl;

			if(header.info_start == 0){

				cerr << "** Info start is 0 -- database is not complete! **" << endl;
				return EXIT_SUCCESS;
			}

			BitVector slice(header.num_filter);

			const size_t num_slice = min(num_bit_slice, header.filter_len() );

			if(num_slice > 0){

				out << "Raw bits for the first " << num_slice << " bitslices" << endl;

				// Read the first few slices from the file
				for(size_t i = 0;i < num_slice;++i){

					slice.read(fin);

					if(!fin){
						throw __FILE__":search: Error reading slice from file (1)";
					}

					out << i;

					for(size_t j = 0;j < header.num_filter;++j){
						out << ' ' << (slice.get_bit(j) ? 1 : 0);
					}

					out << endl;
				}
			}

			// Read and display any annotation information associated with this file
			unsigned long int info_start = header.info_start;

			// Read and throw away the location of each filter annotation
			//unsigned long int *annotation_loc = new unsigned long int [header.num_filter];

			//if(annotation_loc == NULL){
			//	throw __FILE__ ":main: Unable to allocate annotation_loc";
			//}

			//fin.seekg(info_start);
			fin.seekg( info_start + header.num_filter*sizeof(unsigned long int) );

			if(!fin){
				throw __FILE__ ":main: Unable to jump to the start of FilterInfo data";
			}
			//fin.read( (char*)annotation_loc, header.num_filter*sizeof(unsigned long int) );

			//delete [] annotation_loc;

			for(unsigned int i = 0;i < header.num_filter;++i){

				FilterInfo info;

				binary_read(fin, info);

				out << "Annotation information for Bloom filter " << i << endl;

				out << "\trun_accession = " << ( (info.run_accession == INVALID_ACCESSION) ? 
					"NA" : accession_to_str(info.run_accession) ) << endl;
				out << "\texperiment_accession = " << ( (info.experiment_accession == INVALID_ACCESSION) ? 
					"NA" : accession_to_str(info.experiment_accession) ) << endl;
				out << "\texperiment_title = " << (info.experiment_title.empty() ? 
					"NA" : info.experiment_title) << endl;
				out << "\texperiment_design_description = " << (info.experiment_design_description.empty() ? 
					"NA" : info.experiment_design_description) << endl;
				out << "\texperiment_library_name = " << (info.experiment_library_name.empty() ? 
					"NA" : info.experiment_library_name) << endl;
				out << "\texperiment_library_strategy = " << (info.experiment_library_strategy.empty() ? 
					"NA" : info.experiment_library_strategy) << endl;
				out << "\texperiment_library_source = " << (info.experiment_library_source.empty() ? 
					"NA" : info.experiment_library_source) << endl;
				out << "\texperiment_library_selection = " << (info.experiment_library_selection.empty() ? 
					"NA" : info.experiment_library_selection) << endl;
				out << "\texperiment_instrument_model = " << (info.experiment_instrument_model.empty() ?
					"NA" : info.experiment_instrument_model) << endl;
				out << "\tsample_accession = " << ( (info.sample_accession == INVALID_ACCESSION) ? 
					"NA" : accession_to_str(info.sample_accession) ) << endl;
				out << "\tsample_taxa = " << (info.sample_taxa.empty() ? 
					"NA" : info.sample_taxa) << endl;

				if( !info.sample_attributes.empty() ){

					out << "\tsample_attributes:" << endl;

					// Enforce a consistent order on the sample attributes (note that the keys() function
					// sorts the keys in ascending order)
					const vector<string> attrib_key = keys(info.sample_attributes);

					for(vector<string>::const_iterator j = attrib_key.begin();j != attrib_key.end();++j){
						
						MAP<std::string, std::string>::const_iterator attrib_iter = info.sample_attributes.find(*j);

						if( attrib_iter == info.sample_attributes.end() ){
							throw __FILE__ ":main: Unable to lookup sample attribute by key name";
						}							
						
						out << "\t\t" << attrib_iter->first << " = " << attrib_iter->second << endl;
					}
				}

				out << "\tstudy_accession = " << ( (info.study_accession == INVALID_ACCESSION) ? "NA" : 
					accession_to_str(info.study_accession) ) << endl;
				out << "\tstudy_title = " << (info.study_title.empty() ? "NA" : info.study_title) << endl;
				out << "\tstudy_abstract = " << (info.study_abstract.empty() ? "NA" : info.study_abstract) << endl;

				// Separate annotation fields with a blank line
				out << endl;
			}
		}
	}
	catch(const char *error){
		cerr << "Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(const string error){
		cerr << "Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(...){
		cerr << "Caught an unhandled error" << endl;
		return EXIT_FAILURE;
	}
	
	
	return EXIT_SUCCESS;
}
