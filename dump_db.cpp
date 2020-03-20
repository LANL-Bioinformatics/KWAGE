// Dump the contents of BIGSI++ database file in a human readable
// format.
// J. D. Gans
// Bioscience Division, B-11
// Los Alamos National Laboratory
// Mon Mar  2 12:02:07 2020

#include <iostream>
#include <iomanip>
#include <fstream>

#include <stdlib.h>
#include <math.h>

#include "bigsi++.h"
#include "bloom.h"
#include "slice_z.h"

using namespace std;

int main(int argc, char *argv[])
{
	try{
		
		const size_t num_slice_dump = 25;
		
		if(argc != 2){
		
			cerr << "Usage: " << argv[0] 
				<< " <BIGSI++ database file>" << endl;
			return EXIT_SUCCESS;
		}
		
		const string filename = argv[1];
		
		ifstream fin(filename.c_str(), ios::binary);
		
		if(!fin){
		
			cerr << "Unable to open " << filename << " for reading" << endl;
			return EXIT_FAILURE;
		}
		
		DBFileHeader header;
		
		binary_read(fin, header);
		
		if(!fin){
		
			cerr << "Unable to read database header" << endl;
			return EXIT_FAILURE;
		}
		
		cout << "Header information for " << filename << endl;
		cout << "\tmagic = " << header.magic << endl;
		cout << "\tversion = " << header.version << endl;
		cout << "\tcrc32 = " << header.crc32 << endl;
		cout << "\tkmer_len = " << header.kmer_len << endl;
		cout << "\tnum_hash = " << header.num_hash << endl;
		
		const size_t filter_len = header.filter_len();
		
		cout << "\tfilter_len = " << filter_len << endl;
		
		cout << "\tlog_2_filter_len = " << header.log_2_filter_len << endl;
		cout << "\tnum_filter = " << header.num_filter << endl;
		
		switch(header.hash_func){
			case MURMUR_HASH:
				cout << "\thash_func = Murmur" << endl;
				break;
			case UNKNOWN_HASH:
				cout << "\thash_func = Unknown" << endl;
				break;
			default:
				cout << "\thash_func = Invalid" << endl;
				break;
		};
		
		switch(header.compression){
			case NO_COMPRESSION:
				cout << "\tcompression = None" << endl;
				break;
			case ZLIB_RLE_COMPRESSION:
				cout << "\tcompression = zlib (RLE)" << endl;
				break;
			default:
				cout << "\tcompression = Invalid" << endl;
				break;
		};
		
		unsigned char *slice_len = NULL;
		size_t total_compressed_bloom_bytes = 0;

		if(header.compression != NO_COMPRESSION){

			// Read the lengths of the compressed slices
			slice_len = new unsigned char[filter_len];

			if(slice_len == NULL){
				throw __FILE__ ":main: Unable to allocate slice length buffer";
			}

			fin.read( (char*)slice_len, filter_len*sizeof(unsigned char) );

			if(!fin){
				throw __FILE__ ":main: Unable to read slice length buffer";
			}

			for(size_t i = 0;i < filter_len;++i){
				total_compressed_bloom_bytes += slice_len[i];
			}
		}

		// Record the point in the file where the Bloom filter data starts
		const unsigned long int bloom_start = fin.tellg();

		const unsigned int bytes_per_bitslice = header.num_filter/BitVector::BITS_PER_BLOCK +
			( (header.num_filter%BitVector::BITS_PER_BLOCK) > 0 ? 1 : 0);
		
		cout << "There are " << bytes_per_bitslice << " uncompressed bytes per slice" << endl;
		
		InflateSlice<MAX_COMPRESSED_BYTES> decompressor;
		CompressSlice<MAX_COMPRESSED_BYTES> compressor;
		
		BitVector slice(header.num_filter);
		
		float ave_compression = 0.0;
		float stdev_compression = 0.0;

		if(num_slice_dump > 0){

			cout << "Raw bits for the first " << num_slice_dump << " bitslices" << endl;
			cout << setprecision(3);

			// Read the first few slices from the file
			for(size_t i = 0;i < num_slice_dump;++i){

				if( (header.compression == NO_COMPRESSION) || (slice_len[i] == 0) ){

					slice.read(fin);

					if(!fin){
						throw __FILE__":search: Error reading slice from file (1)";
					}
				}
				else{

					// We need to read and inflate the requested bit slice
					decompressor.inflate(fin, slice_len[i]);

					slice.read( decompressor.ptr() );
				}

				// Compress this slice to measure the compression efficiency
				compressor.compress(slice.ptr(), bytes_per_bitslice);
				
				const float fraction_compressed = float(compressor.size())/bytes_per_bitslice;
				
				ave_compression += fraction_compressed;
				stdev_compression += fraction_compressed*fraction_compressed;
				
				cout << i << '\t'<< fraction_compressed << '\t';
				
				for(size_t j = 0;j < header.num_filter;++j){
					cout << slice.get_bit(j);
				}
				
				cout << endl;
			}

			ave_compression /= num_slice_dump;
			stdev_compression = sqrt(stdev_compression/num_slice_dump - ave_compression*ave_compression);
		
			cout << "Average compression: " << ave_compression << " +/- " << stdev_compression << endl;
		}

		if(slice_len != NULL){
			delete [] slice_len;
		}

		// Read and display any annotation information associated with this file
		unsigned long int info_start = 0;

		if(header.compression != NO_COMPRESSION){

			// This is a compressed database file
			info_start = bloom_start + total_compressed_bloom_bytes;
		}
		else{

			// Uncompressed database file
			info_start = bloom_start + filter_len*bytes_per_bitslice;
		}

		// Read and throw away the location of each filter annotation
		unsigned long int *annotation_loc = new unsigned long int [header.num_filter];

		if(annotation_loc == NULL){
			throw __FILE__ ":main: Unable to allocate annotation_loc";
		}

		fin.seekg(info_start);

		fin.read( (char*)annotation_loc, header.num_filter*sizeof(unsigned long int) );

		delete [] annotation_loc;

		for(unsigned int i = 0;i < header.num_filter;++i){

			FilterInfo info;

			binary_read(fin, info);

			cout << "Annotation information for Bloom filter " << i << endl;

			cout << "\trun_accession = " << (info.run_accession.empty() ? 
				"NA" : info.run_accession) << endl;
			cout << "\texperiment_accession = " << (info.experiment_accession.empty() ? 
				"NA" : info.experiment_accession) << endl;
			cout << "\texperiment_title = " << (info.experiment_title.empty() ? 
				"NA" : info.experiment_title) << endl;
			cout << "\texperiment_design_description = " << (info.experiment_design_description.empty() ? 
				"NA" : info.experiment_design_description) << endl;
			cout << "\texperiment_library_name = " << (info.experiment_library_name.empty() ? 
				"NA" : info.experiment_library_name) << endl;
			cout << "\texperiment_library_strategy = " << (info.experiment_library_strategy.empty() ? 
				"NA" : info.experiment_library_strategy) << endl;
			cout << "\texperiment_library_source = " << (info.experiment_library_source.empty() ? 
				"NA" : info.experiment_library_source) << endl;
			cout << "\texperiment_library_selection = " << (info.experiment_library_selection.empty() ? 
				"NA" : info.experiment_library_selection) << endl;
			cout << "\texperiment_instrument_model = " << (info.experiment_instrument_model.empty() ?
				"NA" : info.experiment_instrument_model) << endl;
			cout << "\tsample_accession = " << (info.sample_accession.empty() ? 
				"NA" : info.sample_accession) << endl;
			cout << "\tsample_taxa = " << (info.sample_taxa.empty() ? 
				"NA" : info.sample_taxa) << endl;

			if( !info.sample_attributes.empty() ){

				cout << "\tsample_attributes:" << endl;

				for(MULTIMAP<std::string, std::string>::const_iterator j = info.sample_attributes.begin();
					j != info.sample_attributes.end();++j){
						
						cout << "\t\t" << j->first << " = " << j->second << endl;
				}
			}

			cout << "\tstudy_accession = " << (info.study_accession.empty() ? "NA" : info.study_accession) << endl;
			cout << "\tstudy_title = " << (info.study_title.empty() ? "NA" : info.study_title) << endl;
			cout << "\tstudy_abstract = " << (info.study_abstract.empty() ? "NA" : info.study_abstract) << endl;

			// Separate annotation fields with a blank line
			cout << endl;
		}
	}
	catch(const char *error){
		cerr << "Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(const string error){
		cerr << "Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(...){
		cerr << "Caught an unhandled error" << endl;
		return EXIT_FAILURE;
	}
	
	
	return EXIT_SUCCESS;
}
