#include <iostream>

using namespace std;

struct foo;

void bar(const foo &m_var);

struct foo
{
	int x;

	friend void bar(const foo &m_var);
};

int main(int argc, char* argv[])
{

	try{

		foo A;

		bar(A);
	}
	catch(const char *error){
		cerr << "caught the error: " << error << endl;
		return EXIT_FAILURE;
	}
	catch(...){
		cerr << "caught an unhandled error" << endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

void bar(const foo &m_var)
{

	cerr << "m_var.x = " << m_var.x << endl;
}