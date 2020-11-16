#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <deque>
#include <unordered_map>
#include <algorithm>

#include <stdlib.h>
#include <getopt.h>
#include <signal.h>
#include <math.h>

#include "maestro.h"

using namespace std;

// Global variables for MPI
int mpi_numtasks;
int mpi_rank;

double start_time;

void terminate_program(int m_sig);
string report_run_time();

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
			
	signal( SIGINT, terminate_program );
	signal( SIGTERM, terminate_program );
	signal( SIGSEGV, terminate_program );
			
	start_time = MPI_Wtime();
		
	try{
		
		MaestroOptions opt;
		
		// Only rank 0 reads the command line arguments
		if(mpi_rank == 0){
			opt.load(argc, argv);
		}
		
		broadcast(opt, mpi_rank, 0);
		
		if(opt.quit){
			
			MPI_Finalize();
			return EXIT_SUCCESS;
		}
		
		if(mpi_numtasks < 2){

			cerr << "Please run " << argv[0] << " with more than 1 MPI rank" << endl;

			MPI_Finalize();
			return EXIT_SUCCESS;
		}

		if(mpi_rank == 0){

			maestro_main(opt);

			cerr << report_run_time() << endl;
		}
		else{
			worker_main(opt);
		}
	}
	catch(const char *error){

        cerr << "[" << mpi_rank << "] Caught the error: " << error << endl;
		
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return EXIT_FAILURE;
	}
	catch(const std::exception &error){

        cerr << "[" << mpi_rank << "] Caught the error: " << error.what() << endl;
		
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return EXIT_FAILURE;
	}
	catch(...){

		cerr << "[" << mpi_rank << "] Caught an unhandled error" << endl;

		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		return EXIT_FAILURE;
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}

// Report our rank, the signal we caught and the time spent running
void terminate_program(int m_sig)
{
        cerr << "[" << mpi_rank << "] caught signal " << m_sig << endl;
        cerr << report_run_time() << endl;
        
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

// Run time computes the total run time. The results are formatted as a string.
string report_run_time()
{
        double elapsed_time = MPI_Wtime() - start_time; // In sec
        
        const double elapsed_sec = fmod(elapsed_time, 60.0);
        
        elapsed_time = (elapsed_time - elapsed_sec)/60.0; // In min
        
        const double elapsed_min = fmod(elapsed_time, 60.0);
        elapsed_time = (elapsed_time - elapsed_min)/60.0; // In hour
        
        const double elapsed_hour = fmod(elapsed_time, 24.0);
        elapsed_time = (elapsed_time - elapsed_hour)/24.0; // In day
        
        stringstream sout;
        
        sout << "Run time is " 
                << elapsed_time 
                << " days, "
                << elapsed_hour
                << " hours, "
                << elapsed_min
                << " min and "
                << elapsed_sec
                << " sec";
        
        return sout.str();
}
