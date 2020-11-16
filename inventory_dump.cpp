#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include "bloom.h"
#include "binary_io.h"

using namespace std;

int main(int argc, char *argv[])
{
	try{

		if(argc != 2){

			cerr << "Usage: " << argv[0] << " <binary metadata file>" << endl;
			return EXIT_SUCCESS;
		}

		ifstream fin(argv[1], ios::binary);

		if(!fin){

			cerr << "Unable to open metadata file: " << argv[1] << endl;
			return EXIT_FAILURE;
		}

		FilterInfo info;
		size_t num_info;

		binary_read(fin, num_info);

		cout << "Metadata file contains " << num_info << " FilterInfo objects" << endl;

		for(size_t i = 0;i < num_info;++i){

			binary_read(fin, info);

			if(info.run_accession == INVALID_ACCESSION){
				cout << "Invalid run accession" << endl;
			}
			else{
				cout << accession_to_str(info.run_accession) << endl;
			}

			cout << "\tspots : " << info.number_of_spots << endl;
			cout << "\tbases : " << info.number_of_bases << endl;
			cout << "\tdate_received : " << info.date_received << endl;

			if(info.experiment_accession == INVALID_ACCESSION){
				cout << "\texperiment_accession : Invalid" << endl;
			}
			else{
				cout << "\texperiment_accession : " << accession_to_str(info.experiment_accession) << endl;
			}

			cout << "\texperiment_title : " << info.experiment_title << endl;
			cout << "\texperiment_design_description : " << info.experiment_design_description << endl;
			cout << "\texperiment_library_name : " << info.experiment_library_name << endl;
			cout << "\texperiment_library_strategy : " << info.experiment_library_strategy << endl;
			cout << "\texperiment_library_source : " << info.experiment_library_source << endl;
			cout << "\texperiment_library_selection : " << info.experiment_library_selection << endl;
			cout << "\texperiment_instrument_model : " << info.experiment_instrument_model << endl;

			if(info.sample_accession == INVALID_ACCESSION){
				cout << "\tsample_accession : Invalid" << endl;
			}
			else{
				cout << "\tsample_accession : " << accession_to_str(info.sample_accession) << endl;
			}

			cout << "\tsample_taxa : " << info.sample_taxa << endl;

			if( !info.sample_attributes.empty() ){

				cout << "\tsample_attributes :" << endl;

				for(MAP<string, string>::const_iterator j = info.sample_attributes.begin();
					j != info.sample_attributes.end();++j){
					
					cout << "\t\t" << j->first << " : " << j->second << endl;
				}
			}

			if(info.study_accession == INVALID_ACCESSION){
				cout << "\tstudy_accession : Invalid" << endl;
			}
			else{
				cout << "\tstudy_accession : " << accession_to_str(info.study_accession) << endl;
			}
			
			cout << "\tstudy_title : " << info.study_title << endl;
			cout << "\tstudy_abstract : " << info.study_abstract << endl;
		}
	}
	catch(const char *error){

		cerr << "Caught the error: " << error << endl;
		return EXIT_FAILURE;
	}
	catch(...){

		cerr << "Caught an unhandled error!" << endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
