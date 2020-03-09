#include "update.h"
#include <iostream>

using namespace std;

void UpdateInfo::flush()
{
	if(valid){
	
		// Clear the buffer
		for(size_t i = 0;i < buffer_size;i++){
			cerr << '\b';
		}

		for(size_t i = 0;i < buffer_size;i++){
			cerr << ' ';
		}

		for(size_t i = 0;i < buffer_size;i++){
			cerr << '\b';
		}

		const string tmp = str();

		buffer_size = tmp.size();

		cerr << tmp;

		if( flog.is_open() ){
			flog << tmp << endl;
		}

		// Clear the stringstream buffer
		str( string() );
	}
}

void UpdateInfo::close()
{
	if(valid){
	
		cerr << endl;

		// Clear the stringstream buffer
		str( string() );
	}
}
