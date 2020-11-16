#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{

	for(unsigned char x = 0;x < 8;++x){
		cerr << "x = " << int(x) << " --> " << int( ~x & 7 ) << endl;
	}

	return EXIT_SUCCESS;
}