#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>
#include <ncbi-vdb/NGS.hpp> // For openReadCollection
#include <getopt.h>

using namespace std;

#define		SRACAT_VERSION		"0.1"

void clear_buffer(ostream &m_out, const string &m_buffer);

int main(int argc, char *argv[])
{
	try{

		const char* options = "o:v?h";
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
		string output_filename;
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
					output_filename = optarg;
					break;
				case 'h':
				case '?':

					print_usage = true;
					break;
				case 'v':
					verbose = true;
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
			cerr << "\t[-o <output file>] (default is stdout)" << endl;
			cerr << "\t[--qual] (fastq output)" << endl;
			cerr << "\t<SRA accession/file 1> ..." << endl;

			return EXIT_SUCCESS;
		}
		
		ofstream fout;

		if( !output_filename.empty() ){

			fout.open( output_filename.c_str() );

			if(!fout){

				cerr << "Unable to open " << output_filename << " for writing" << endl;
				return EXIT_FAILURE;
			}
		}

		ostream &out = fout.is_open() ? fout : cout;

		for(deque<string>::const_iterator acc = sra_records.begin();acc != sra_records.end();++acc){

			ngs::ReadCollection run( ncbi::NGS::openReadCollection(*acc) );
			
			size_t output_read_count = 0; // For labeling the output reads

			// Following the advice of Kurt Rodamer @ NCBI
			// ** Note that this approach will *miss* a small number of reads that are only
			// ** partially aligned (i.e. only one read of a pair aligned). If all reads
			// ** in an SRA record are needed, then iterate through all reads using
			// ** ngs::ReadIterator( run.getReadRange ( 1, num_read, ngs::Read::all ) );
			// Step 1: Does the SRA record contain aligned reads?
			const size_t num_primary_align = run.getAlignmentCount(ngs::Alignment::primaryAlignment);

			const string accession_name = run.getName();

			if(verbose){
				cerr << accession_name << ": ";
			}

			string verbose_buffer;

			if(num_primary_align > 0){

				if(qual){
					cerr << "fastq output is not supported for aligned reads, switching to fasta" << endl;
				}

				// Read the primaryAlignment
				ngs::AlignmentIterator align_iter = run.getAlignments(ngs::Alignment::primaryAlignment);

				size_t align_count = 0;

				const size_t update_every = max(num_primary_align/100ULL, 1ULL);

				while(align_iter.nextAlignment() ){

					++output_read_count;
					++align_count;

					const ngs::StringRef &seq = align_iter.getAlignedFragmentBases();

					// Only FASTA is supported for now
					out << ">" << accession_name << "." << output_read_count << '\n' << seq << endl;

					if( verbose && (align_count%update_every == 0) ){

						stringstream ssout;

						ssout << "primary " << (100.0*align_count)/num_primary_align << "%";

						clear_buffer(cerr, verbose_buffer);

						verbose_buffer = ssout.str();

						cerr << verbose_buffer;
					}
				}

				// Read the unaligned sequences
				const size_t num_unaligned_read = run.getReadCount(ngs::Read::unaligned);

				if(num_unaligned_read > 0){

					// Need to use getReads() -- getReadRange() does not appear to work for unaligned reads
					ngs::ReadIterator run_iter = ngs::ReadIterator( run.getReads(ngs::Read::unaligned) );

					size_t read_count = 0;

					const size_t update_every = max(num_unaligned_read/100ULL, 1ULL);

					while( run_iter.nextRead() ){
						
						++output_read_count;
						++read_count;

						size_t seq_count = 0;

						while( run_iter.nextFragment() ){
							
							++seq_count;
													
							const ngs::StringRef &seq = run_iter.getFragmentBases();
							
							// Only fasta is supported for now
							out << ">" << accession_name << '.' << output_read_count << '.' << seq_count << '\n' << seq << endl;
						}

						if( verbose && (read_count%update_every == 0) ){

							stringstream ssout;

							ssout << "unaligned " << (100.0*read_count)/num_unaligned_read << "%";

							clear_buffer(cerr, verbose_buffer);

							verbose_buffer = ssout.str();

							cerr << verbose_buffer;
						}
					}				
				}
			}
			else{

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
						
						if(qual){

							// fastq output
							const ngs::StringRef &phred_q = run_iter.getFragmentQualities();

							out << "@" << accession_name << '.' << output_read_count << '.' << seq_count << '\n' << seq << '\n';
							out << '+' << '\n';
							out << phred_q << endl;
						}
						else{

							// fasta output
							out << ">" << accession_name << '.' << output_read_count << '.' << seq_count << '\n' << seq << endl;
						}
					}

					if( verbose && (read_count%update_every == 0) ){

						stringstream ssout;

						ssout << (100.0*read_count)/num_read << "%";

						clear_buffer(cerr, verbose_buffer);

						verbose_buffer = ssout.str();

						cerr << verbose_buffer;
					}
				}
			}

			if(verbose){
				cerr << endl;
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