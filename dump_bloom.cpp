// Dump the contents of CALDERA Bloom filter file in a human readable
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

#include "caldera.h"
#include "bloom.h"

using namespace std;

int main(int argc, char *argv[])
{
	try{
		
		if(argc != 2){
		
			cerr << "Usage: " << argv[0] 
				<< " <CALDERA Bloom filter file>" << endl;
			return EXIT_SUCCESS;
		}
		
		const string filename = argv[1];
		
		ifstream fin(filename.c_str(), ios::binary);
		
		if(!fin){
		
			cerr << "Unable to open " << filename << " for reading" << endl;
			return EXIT_FAILURE;
		}
		
		BloomFilter filter;
		
		binary_read(fin, filter);
		
		if(!fin){
		
			cerr << "Unable to read Bloom filter header" << endl;
			return EXIT_FAILURE;
		}
		
		const BloomParam &param = filter.get_param();

		const size_t filter_len = param.filter_len();

		cout << "Header information for " << filename << endl;
		cout << "\tcrc32 = " << std::hex << filter.get_crc32() << std::dec << endl;
		cout << "\tlength = " << param.filter_len() << endl;
		cout << "\tlog_2 length = " << param.log_2_filter_len << endl;
		cout << "\tnum_hash = " << param.num_hash << endl;
		cout << "\tkmer_len = " << param.kmer_len << endl;
				
		switch(param.hash_func){
			case MURMUR_HASH_32:
				cout << "\thash_func = Murmur32" << endl;
				break;
			case UNKNOWN_HASH:
				cout << "\thash_func = Unknown" << endl;
				break;
			default:
				cout << "\thash_func = Invalid" << endl;
				break;
		};

		const FilterInfo &info = filter.get_info();

		cout << "Annotation information for Bloom filter " << endl;

		cout << "\trun_accession = " << ( (info.run_accession == INVALID_ACCESSION) ? 
			"NA" : accession_to_str(info.run_accession) ) << endl;
		cout << "\texperiment_accession = " << ( (info.experiment_accession == INVALID_ACCESSION) ? 
			"NA" : accession_to_str(info.experiment_accession) ) << endl;
		cout << "\texperiment_title = " << (info.experiment_title.empty() ? 
			"NA" : info.experiment_title) << endl;
		cout << "\texperiment_design_description = " << (info.experiment_design_description.empty() ? 
			"NA" : info.experiment_design_description) << endl;
		cout << "\texperiment_library_name = " << (info.experiment_library_name.empty() ? 
			"NA" : info.experiment_library_name) << endl;
		cout << "\texperiment_library_strategy = " << (info.experiment_library_strategy.empty() ? 
			"NA" : info.experiment_library_strategy) << endl;
		cout << "\texperiment_library_source = " << (info.experiment_library_source.empty() ? 
			"NA" : info.experiment_library_source) << endl;
		cout << "\texperiment_library_selection = " << (info.experiment_library_selection.empty() ? 
			"NA" : info.experiment_library_selection) << endl;
		cout << "\texperiment_instrument_model = " << (info.experiment_instrument_model.empty() ?
			"NA" : info.experiment_instrument_model) << endl;
		cout << "\tsample_accession = " << ( (info.sample_accession == INVALID_ACCESSION) ? 
			"NA" : accession_to_str(info.sample_accession) ) << endl;
		cout << "\tsample_taxa = " << (info.sample_taxa.empty() ? 
			"NA" : info.sample_taxa) << endl;

		if( !info.sample_attributes.empty() ){

			cout << "\tsample_attributes:" << endl;

			for(MAP<std::string, std::string>::const_iterator j = info.sample_attributes.begin();
				j != info.sample_attributes.end();++j){
					
					cout << "\t\t" << j->first << " = " << j->second << endl;
			}
		}

		cout << "\tstudy_accession = " << ( (info.study_accession == INVALID_ACCESSION) ? "NA" : 
			accession_to_str(info.study_accession) ) << endl;
		cout << "\tstudy_title = " << (info.study_title.empty() ? "NA" : info.study_title) << endl;
		cout << "\tstudy_abstract = " << (info.study_abstract.empty() ? "NA" : info.study_abstract) << endl;
		
		cout << "Raw bits:" << endl;

		for(size_t i = 0;i < filter_len;++i){
			cout << i << '\t' << (filter.get_bit(i) ? 1 : 0) << endl;
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
