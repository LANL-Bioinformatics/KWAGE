#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include "bloom.h"
#include "binary_io.h"

using namespace std;

void parse_accessions(ifstream &m_fin, vector<SraAccession> &m_acc);

int main(int argc, char *argv[])
{
	try{
		if(argc != 3){

			cerr << "Usage: " << argv[0] << " <binary metadata file 1> <binary metadata file 2>" << endl;
			return EXIT_SUCCESS;
		}

		ifstream fin_1(argv[1], ios::binary);

		if(!fin_1){

			cerr << "Unable to open file 1: " << argv[1] << endl;
			return EXIT_FAILURE;
		}

		ifstream fin_2(argv[2], ios::binary);

		if(!fin_2){

			cerr << "Unable to open file 2: " << argv[2] << endl;
			return EXIT_FAILURE;
		}

		vector<SraAccession> acc_1;
		vector<SraAccession> acc_2;

		cerr << "Reading file 1: " << argv[1] << endl;

		try{
			parse_accessions(fin_1, acc_1);
			sort( acc_1.begin(), acc_1.end() );
		}
		catch(...){
			cerr << "Unable to parse file 1: " << argv[1] << endl;
		}

		cerr << "Reading file 2: " << argv[2] << endl;

		try{

			parse_accessions(fin_2, acc_2);
			sort( acc_2.begin(), acc_2.end() );
		}
		catch(...){
			cerr << "Unable to parse file 2: " << argv[2] << endl;
		}

		vector<SraAccession>::const_iterator iter_1 = acc_1.begin();
		vector<SraAccession>::const_iterator iter_2 = acc_2.begin();

		cerr << "Comparing accession sets" << endl;

		while(true){

			if( iter_1 == acc_1.end() ){

				cout << "Reached the last accession of the first file" << endl;
				cout << "There are " << (acc_2.size() - ( iter_2 - acc_2.begin() ) ) 
					<< " accessions remaining in the second file" << endl;
				break;
			}

			if( iter_2 == acc_2.end() ){

				cout << "Reached the last accession of the second file" << endl;
				cout << "There are " << (acc_1.size() - ( iter_1 - acc_1.begin() ) ) 
					<< " accessions remaining in the first file" << endl;
				break;
			}

			if(*iter_1 < *iter_2){

				cout << "1: " << accession_to_str(*iter_1) << endl;

				++iter_1;
				continue;
			}

			if(*iter_2 < *iter_1){
				
				cout << "2: " << accession_to_str(*iter_2) << endl;
				
				++iter_2;
				continue;
			}

			// *iter_1 == *iter_2
			++iter_1;
			++iter_2;
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

void parse_accessions(ifstream &m_fin, vector<SraAccession> &m_acc)
{
	FilterInfo info;
	size_t num_info;

	binary_read(m_fin, num_info);

	m_acc.resize(num_info);

	for(size_t i = 0;i < num_info;++i){

		binary_read(m_fin, info);

		m_acc[i] = info.run_accession;
	}
}