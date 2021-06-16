#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <deque>
#include <ncbi-vdb/NGS.hpp> // For openReadCollection
#include <getopt.h>
#include <zlib.h>
#include <stdio.h>

using namespace std;

#define		SRACAT_VERSION		"0.2"

void clear_buffer(ostream &m_out, const string &m_buffer);

int main(int argc, char *argv[])
{
	try{

		const char* options = "o:zv?h";
		int config_opt = 0;
		int long_index = 0;

		struct option long_opts[] = {
			{"qual", false, &config_opt, 1},
			{0,0,0,0} // Terminate options list
		};

		int opt_code;
		opterr = 0;

		bool print_usage = false;
		bool verbose = false;
		bool qual = false; // fasta by default, but switch to fastq if quality scores are requested
		bool compress_output = false;
		string output_prefix;
		deque<string> sra_records;

		while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

			switch( opt_code ){
				case 0:

					if(config_opt == 1){ // qual
				
						qual = true;
						break;
					}
					
					cerr << "Unknown command line flag!" << endl;
					return EXIT_FAILURE;
				case 'o':
					output_prefix = optarg;
					break;
				case 'h':
				case '?':

					print_usage = true;
					break;
				case 'v':
					verbose = true;
					break;
				case 'z':
					compress_output = true;
					break;
				default:
					cerr << '\"' << (char)opt_code << "\" is not a valid option!" << endl;
					return EXIT_FAILURE;
			};
		}

		for(int i = optind;i < argc;++i){
			sra_records.push_back(argv[i]);
		}

		if(print_usage){

			cerr << "Usage for sracat (v. " << SRACAT_VERSION << ")" << endl;
			cerr << "\t[-o <output file *prefix*>] (default is stdout)" << endl;
			cerr << "\t[-z] (zlib-based compression of file-based output; default is no compression)" << endl;
			cerr << "\t[--qual] (fastq output)" << endl;
			cerr << "\t<SRA accession/file 1> ..." << endl;

			return EXIT_SUCCESS;
		}
		
		// When writing to a file, each "spot" can be associated with multiple reads (i.e. read 1 and 2 for Illumina)
		// Partition the different read types into separate files
		deque<gzFile> fout;

		for(deque<string>::const_iterator acc = sra_records.begin();acc != sra_records.end();++acc){

			ngs::ReadCollection run( ncbi::NGS::openReadCollection(*acc) );
			
			size_t output_read_count = 0; // For labeling the output reads

			const string accession_name = run.getName();

			if(verbose){
				cerr << accession_name << ": ";
			}

			string verbose_buffer;

			const size_t num_read = run.getReadCount(ngs::Read::all);

			ngs::ReadIterator run_iter = 
				ngs::ReadIterator( run.getReadRange ( 1, num_read, ngs::Read::all ) );

			size_t read_count = 0;

			const size_t update_every = max(num_read/100ULL, 1ULL);

			while( run_iter.nextRead() ){
				
				++output_read_count;
				++read_count;

				size_t seq_count = 0;

				while( run_iter.nextFragment() ){
					
					++seq_count;

					const ngs::StringRef &seq = run_iter.getFragmentBases();
					
					gzFile out = NULL;

					if( !output_prefix.empty() ){

						if(fout.size() < seq_count){

							// Open a new file
							stringstream filename;
							
							filename << output_prefix << '_' << seq_count;

							if(qual){
								filename << ".fastq";
							}
							else{
								filename << ".fna";
							}

							if(compress_output){

								filename << ".gz";
								fout.push_back( gzopen( filename.str().c_str(), "w") );
							}
							else{
								// THe 'T' flag disables compression
								fout.push_back( gzopen( filename.str().c_str(), "wT") );
							}

							if(fout.back() == NULL){
								
								cerr << "Unable to open: " << filename.str() << endl;
								throw "Unable to open output file for writing";
							}
						}

						out = fout[seq_count - 1];
					}

					if(qual){

						// fastq output
						const ngs::StringRef &phred_q = run_iter.getFragmentQualities();

						if(out == NULL){
							printf( "@%s.%lu.%lu\n%s\n+\n%s\n", accession_name.c_str(), output_read_count, seq_count, 
								seq.toString().c_str(), phred_q.toString().c_str() );
						}
						else{
							gzprintf( out, "@%s.%lu\n%s\n+\n%s\n", accession_name.c_str(), output_read_count, 
								seq.toString().c_str(), phred_q.toString().c_str() );
						}

						//out << "@" << accession_name << '.' << output_read_count << '.' << seq_count << '\n' << seq << '\n';
						//out << '+' << '\n';
						//out << phred_q << endl;
					}
					else{

						// fasta output
						if(out == NULL){
							printf( ">%s.%lu.%lu\n%s\n", accession_name.c_str(), output_read_count, seq_count, 
								seq.toString().c_str() );
						}
						else{
							gzprintf( out, ">%s.%lu\n%s\n", accession_name.c_str(), output_read_count, 
								seq.toString().c_str() );
						}

						//out << ">" << accession_name << '.' << output_read_count << '.' << seq_count << '\n' << seq << endl;
					}
				}

				if( verbose && (read_count%update_every == 0) ){

					stringstream ssout;

					ssout << setprecision(3) << (100.0*read_count)/num_read << "%";

					clear_buffer(cerr, verbose_buffer);

					verbose_buffer = ssout.str();

					cerr << verbose_buffer;
				}
			}

			if(verbose){
				cerr << endl;
			}
		}

		for(deque<gzFile>::iterator f = fout.begin();f != fout.end();++f){

			if(*f != NULL){

				gzclose(*f);
				*f = NULL;
			}
		}
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

void clear_buffer(ostream &m_out, const string &m_buffer)
{
	for(string::const_iterator i = m_buffer.begin();i != m_buffer.end();++i){
		m_out << '\b';
	}

	for(string::const_iterator i = m_buffer.begin();i != m_buffer.end();++i){
		m_out << ' ';
	}

	for(string::const_iterator i = m_buffer.begin();i != m_buffer.end();++i){
		m_out << '\b';
	}
}