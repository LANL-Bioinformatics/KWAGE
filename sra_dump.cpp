#include <iostream>

// DEBUG
#include "mem_usage.h"

#include <ncbi-vdb/NGS.hpp> // For openReadCollection

using namespace std;

int main(int argc, char *argv[])
{
	try{
		if(argc != 2){

			cerr << "Usage: " << argv[0] << " <SRA accession>" << endl;
			return EXIT_SUCCESS;
		}

		const string acc = argv[1];

		ngs::ReadCollection run(  ncbi::NGS::openReadCollection(acc) );
		
		// Following the advice of Kurt Rodamer @ NCBI
		// ** Note that this approach will *miss* a small number of reads that are only
		// ** partially aligned (i.e. only one read of a pair aligned). If all reads
		// ** in an SRA record are needed, then iterate through all reads using
		// ** ngs::ReadIterator( run.getReadRange ( 1, num_read, ngs::Read::all ) );
		// Step 1: Does the SRA record contain aligned reads?
		const size_t num_primary_align = run.getAlignmentCount(ngs::Alignment::primaryAlignment);

		if(num_primary_align > 0){

			// Step 2: Read the primaryAlignment
			ngs::AlignmentIterator align_iter = run.getAlignments(ngs::Alignment::primaryAlignment);

			size_t align_count = 0;

			// DEBUG
			const size_t update_every = max(num_primary_align/100ULL, 1ULL);
			time_t profile = time(NULL);

			while(align_iter.nextAlignment() ){

				++align_count;

				const ngs::StringRef &seq = align_iter.getAlignedFragmentBases();

				cout << align_count << '\t' << seq << endl;

				if(align_count%update_every == 0){

					profile = time(NULL) - profile;

					cerr << "primary " << (100.0*align_count)/num_primary_align << "% complete @ " 
						<< double(update_every)/profile<< " alignments/sec; " << memory_usage() << "% memory usage" << endl;
					
					profile = time(NULL);
				}
			}

			// Step 3: Read the unaligned sequences
			const size_t num_unaligned_read = run.getReadCount(ngs::Read::unaligned);

			if(num_unaligned_read > 0){

				// Need to use getReads() -- getReadRange() does not appear to work for unaligned reads
				ngs::ReadIterator run_iter = ngs::ReadIterator( run.getReads(ngs::Read::unaligned) );

				size_t read_count = 0;

				// DEBUG
				const size_t update_every = max(num_unaligned_read/100ULL, 1ULL);
				time_t profile = time(NULL);

				while( run_iter.nextRead() ){
					
					++read_count;

					size_t seq_count = 0;

					while( run_iter.nextFragment() ){
						
						++seq_count;
												
						const ngs::StringRef &seq = run_iter.getFragmentBases();
						
						//cout << read_count << '.' << seq_count << '\t' << seq << endl;
						cout << read_count << '.' << seq_count << '\t' << seq << endl;
					}

					if(read_count%update_every == 0){

						profile = time(NULL) - profile;

						cerr << "unaligned " << (100.0*read_count)/num_unaligned_read << "% complete @ " 
							<< double(update_every)/profile<< " reads/sec; " << memory_usage() << "% memory usage" << endl;
						
						profile = time(NULL);
					}
				}				
			}
		}
		else{

			const size_t num_read = run.getReadCount(ngs::Read::all);

			ngs::ReadIterator run_iter = 
				ngs::ReadIterator( run.getReadRange ( 1, num_read, ngs::Read::all ) );

			size_t read_count = 0;

			// DEBUG
			const size_t update_every = max(num_read/100ULL, 1ULL);
			time_t profile = time(NULL);

			while( run_iter.nextRead() ){
				
				++read_count;

				size_t seq_count = 0;

				while( run_iter.nextFragment() ){
					
					++seq_count;

					const ngs::StringRef &seq = run_iter.getFragmentBases();
					
					cout << read_count << '.' << seq_count << '\t' << seq << endl;
				}

				if(read_count%update_every == 0){

					profile = time(NULL) - profile;

					cerr << "raw " << (100.0*read_count)/num_read << "% complete @ " 
						<< double(update_every)/profile<< " reads/sec; " << memory_usage() << "% memory usage" << endl;
					
					profile = time(NULL);
				}
			}
		}

		#ifdef ORIGINAL_VERSION
		// Note that num_read is the number of either paired or
		// unpaired reads. For paired reads, this is half the
		// the number of sequences!
		const size_t num_read = run.getReadCount(ngs::Read::all);

		// DEBUG
		cerr << "Found " << num_read << " reads in " << acc << endl;

		ngs::ReadIterator run_iter = 
			ngs::ReadIterator( run.getReadRange ( 1, num_read, ngs::Read::all ) );

		size_t read_count = 0;

		// DEBUG
		const size_t update_every = 1000000;
		time_t profile = time(NULL);

		while( run_iter.nextRead() ){
			
			++read_count;

			size_t seq_count = 0;

			while( run_iter.nextFragment() ){
				
				++seq_count;
										
				//const string seq = run_iter.getFragmentBases().toString();
				const ngs::StringRef &seq = run_iter.getFragmentBases();
				
				//cout << read_count << '.' << seq_count << '\t' << seq << endl;
				cout << read_count << '.' << seq_count << '\t' << seq << endl;
			}

			// DEBUG
			if(read_count%update_every == 0){
				
				profile = time(NULL) - profile;

				cerr << (100.0*read_count)/num_read << "% complete @ " 
					<< double(update_every)/profile<< " reads/sec; " << memory_usage() << "% memory usage" << endl;
				
				profile = time(NULL);
			}
		}
		#endif // ORIGINAL_VERSION

	}
	catch(const char* error){

		cerr << "Caught the error: " << error << endl;
		return EXIT_FAILURE;
	}
	catch(...){

		cerr << "Caught an unhandled error" << endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}