#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

int main(int argc, char *argv[])
{
	try{

		const string name = "binary_debug.bin";

		const size_t buffer_size = 2048;

		unsigned char *buffer = new unsigned char[buffer_size];

		if(!buffer){
			throw __FILE__ ":main: Unable to allocate buffer";
		}

		memset(buffer, 0, buffer_size);

		ofstream fout(name.c_str(), ios::binary);

		if(!fout){
			throw __FILE__ ":main: Unable to open outoput file";
		}

		delete [] buiffer;
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
