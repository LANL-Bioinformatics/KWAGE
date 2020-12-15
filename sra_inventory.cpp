// Download SRA sequence and metadata
#include <iostream>
#include <fstream>
#include <algorithm> // for sort()

#include "options.h"
#include "bloom.h"
#include "parse_tar.h"
#include "file_util.h"
#include "kwage.h" // For INVENTORY_VERSION
#include "split.h"
#include "date.h"
#include "string_conversion.h"
#include "binary_io.h"

using namespace std;

typedef enum{
	RUN_XML,
	EXPERIMENT_XML,
	SAMPLE_XML,
	STUDY_XML,
	//SUBMISSION_XML,
	//ANALYSIS_XML,
	SRA_RUN_MEMBERS,
	SRA_ACCESSIONS,
	UNKNOWN
} SRAFileType;

#define		SEC_PER_DAY		(60*60*24)

// The following information is all extracted from the same tab-delimited file
struct AccessionInfo
{
	Date d;
	size_t num_spots;
	size_t num_bases;

	AccessionInfo()
	{
		// Do nothing
	};

	AccessionInfo(const Date &m_date, const size_t &m_num_spots, const size_t &m_num_bases) :
		d(m_date), num_spots(m_num_spots), num_bases(m_num_bases)
	{

	};
};

SRAFileType sra_file_type(const string &m_filename);
string parse_key(const string &m_buffer, string m_key);
string parse_xml(const string &m_key, const string &m_buffer);
void parse_sra_text(deque<FilterInfo> &m_db, const string &m_metadata_file);
void parse_sra_metadata(deque<FilterInfo> &m_db, MAP<SraAccession, char*> &m_sample_attributes,
	const string &m_metadata_file,
	const InventoryOptions &m_opt);
void merge_db(deque<FilterInfo> &m_db, MAP<SraAccession, char*> &m_sample_attributes,
	const MAP<SraAccession, FilterInfo> &m_xml_info, 
	size_t &m_num_experiment, size_t &m_num_sample, size_t &m_num_study);
MAP<string, string> buffer_to_map(char *m_buffer);
char* map_to_buffer(const MAP<string, string> &m_map);

int main(int argc, char* argv[])
{
	try{

		InventoryOptions opt(argc, argv);

		if(opt.quit){
			return EXIT_SUCCESS;
		}

		time_t profile = time(NULL);
		
		// Step 0: Open a file for writing
		ofstream fout;

		if(!opt.list_only){

			fout.open(opt.output_file.c_str(), ios::binary);

			if(!fout){

				cerr << "Unable to open output file for writing binary inventory" << endl;
				return EXIT_FAILURE;
			}
		}

		// Step 1: Parse the SRA metadata to extract all of the SRA runs and selected metadata
		deque<FilterInfo> _db;
		MAP<SraAccession, char*> sample_attributes;

		cerr << "Reading NCBI metadata from " << opt.metadata_file << endl;

		time_t profile_parse_xml = time(NULL);

		parse_sra_metadata(_db, sample_attributes, opt.metadata_file, opt);

		profile_parse_xml = time(NULL) - profile_parse_xml;
			
		size_t num_sra_run = _db.size();

		cout << "Found metadata for " << num_sra_run << " SRA runs in " 
			<< profile_parse_xml << " sec" << endl;

		// Apply the required source filter
		if( !opt.required_source.empty() ){

			cerr << "Applying the requested source filter" << endl;

			size_t num_valid = 0;
			size_t num_invalid = 0;

			for(deque<FilterInfo>::iterator r = _db.begin();r != _db.end();++r){

				// Exclude SRA accessions that do not have the required experimental source
				// This search is *case sensistive* for now
				if( opt.required_source.find(r->experiment_library_source) == opt.required_source.end() ){						

					if(r->valid){

						r->valid = false;
						++num_invalid;
					}
				}
				else{
					if(r->valid){
						++num_valid;
					}
				}
			}

			cerr << "\tFiltered out " << num_invalid << " SRA run accessions" << endl;
			cerr << "\tKept " << num_valid << " SRA run accessions" << endl;
		}

		// Apply the required strategy filter
		if( !opt.required_strategy.empty() ){

			cerr << "Applying the requested strategy filter" << endl;
			size_t num_valid = 0;
			size_t num_invalid = 0;

			for(deque<FilterInfo>::iterator r = _db.begin();r != _db.end();++r){

				// Exclude SRA accessions that do not have the required experimental strategy
				// This search is *case sensistive* for now
				if( opt.required_strategy.find(r->experiment_library_strategy) == opt.required_strategy.end() ){
					
					if(r->valid){

						r->valid = false;
						++num_invalid;
					}
				}
				else{
					if(r->valid){
						++num_valid;
					}
				}
			}

			cerr << "\tFiltered out " << num_invalid << " SRA run accessions" << endl;
			cerr << "\tKept " << num_valid << " SRA run accessions" << endl;
		}

		cerr << "Found sample attributes for " << sample_attributes.size() << " sample accessions" << endl;

		// Apply the data filter
		size_t num_date_invalid = 0;

		for(deque<FilterInfo>::iterator r = _db.begin();r != _db.end();++r){

			if( (r->date_received < opt.begin_date) || (r->date_received > opt.end_date) ){

				if(r->valid){

					r->valid = false;
					++num_date_invalid;
				}
			}
		}

		cerr << "Removed " << num_date_invalid << " SRA run accessions that did not satisfy the date criteria" << endl;

		// Apply the requested inclusion filter
		if( !opt.include_accessions.empty() ){

			cerr << "Applying the requested inclusion filter" << endl;

			size_t num_valid = 0;
			size_t num_invalid = 0;

			for(deque<FilterInfo>::iterator r = _db.begin();r != _db.end();++r){

				// Exclude SRA accessions that are not in the set of accessions to include
				deque<SraAccession>::const_iterator iter = 
					lower_bound(opt.include_accessions.begin(), opt.include_accessions.end(), r->run_accession);
				
				if( (iter == opt.include_accessions.end() ) || (*iter != r->run_accession) ){
					
					if(r->valid){

						r->valid = false;
						++num_invalid;
					}
				}
				else{
					if(r->valid){
						++num_valid;
					}
				}
			}

			cerr << "\tFiltered out " << num_invalid << " SRA run accessions" << endl;
			cerr << "\tKept " << num_valid << " SRA run accessions (out of a requested list of " 
				<< opt.include_accessions.size() << " run accessions)" << endl;
		}

		cerr << "Repacking SRA records" << endl;

		size_t num_zero_base = 0; // This reflects missing data
		size_t total_bases = 0;
		size_t total_spots = 0;

		deque<FilterInfo> db;

		while( !_db.empty() ){

			// Only keep the valid SRA runs
			if(_db.front().valid){

				FilterInfo& ref = _db.front();

				total_bases += ref.number_of_bases;
				total_spots += ref.number_of_spots;

				db.push_back( ref );

				if(ref.number_of_bases == 0){
				
					++num_zero_base;

					// Provide a warning for SRA records that don't have
					// base or spot counts. However, don't remove these records as the
					// absense of base & spot counts is only a formatting issue and most of these
					// are still valid records.
					//cout << accession_to_str(ref.run_accession) << " has zero bases and " 
					//	<< ref.number_of_spots << " spots: " 
					//	<< ref.date_received << endl;
				}
			}

			_db.pop_front();
		}

		num_sra_run = db.size();

		cout << "There are " << num_sra_run << " valid SRA records" << endl;

		// KB = 1024 = 2^10 == 1 << 10
		// MB = 1 << 20
		// GB = 1 << 30
		// TB = 1 << 40
		// PB = 1 << 50
		const double petabases = double(total_bases)/(size_t(1) << 50);

		cout << "Total number of bases = " << total_bases << " = " << petabases << " Petabases" << endl;
		cout << "Total number of spots = " << total_spots << endl;

		cerr << "Sorting SRA records by number of bases ... ";
		sort( db.begin(), db.end() );
		cerr << "done." << endl;


		if( !db.empty() ){

			size_t min_bases = 0;
			SraAccession smallest_sra;

			for(deque<FilterInfo>::const_iterator i = db.begin();i != db.end();++i){
				
				// Find the SRA record with the smallest, non-zero number of bases
				if(i->number_of_bases > 0){
					
					smallest_sra = i->run_accession;
					min_bases = i->number_of_bases;
					break;
				}
			}
	
			const size_t max_bases = db.back().number_of_bases;

			cout << "The number of bases range from " << min_bases << " to " << max_bases << endl;

			cout << "The smallest SRA record is " << accession_to_str(smallest_sra) << endl;
			cout << "The largest SRA record is " << accession_to_str( db.back().run_accession) << endl;

			cout << "The average number of bases/run = " << total_bases/db.size() << endl;
			cout << "The average number of spots/run = " << total_spots/db.size() << endl;

			#ifdef MAKE_HISTOGRAM
			// Make a histogram of bases/record
			const size_t width = (max_bases - min_bases);

			if(width > 0){

				// Use Gigabase bin size for now
				const size_t num_bin = max( size_t(1), (width)/(size_t(1) << 30) );

				vector<size_t> hist(num_bin);

				for(deque<FilterInfo>::const_iterator i = db.begin();i != db.end();++i){

					const size_t index = min( size_t( num_bin*double(i->number_of_bases - min_bases)/width ), 
						num_bin - 1);

					++hist[index];
				}

				cout << "Histogram of base counts" << endl;

				for(size_t i = 0;i < num_bin;++i){
					
					if(hist[i] > 0){
						cout << (i + 1)*width/num_bin << '\t' << hist[i] << endl;
					}
				}
			}
			#endif // MAKE_HISTOGRAM
		}

		if(opt.list_only){

			while( !db.empty() ){
				
				FilterInfo i = db.front();

				db.pop_front();

				cout << accession_to_str(i.run_accession) << endl;
				cout << "\tspots : " << i.number_of_spots << endl;
				cout << "\tbases : " << i.number_of_bases << endl;
				cout << "\tdate_received : " << i.date_received << endl;
				cout << "\texperiment_accession : " << accession_to_str(i.experiment_accession) << endl;
				cout << "\texperiment_title : " << i.experiment_title << endl;
				cout << "\texperiment_design_description : " << i.experiment_design_description << endl;
				cout << "\texperiment_library_name : " << i.experiment_library_name << endl;
				cout << "\texperiment_library_strategy : " << i.experiment_library_strategy << endl;
				cout << "\texperiment_library_source : " << i.experiment_library_source << endl;
				cout << "\texperiment_library_selection : " << i.experiment_library_selection << endl;
				cout << "\texperiment_instrument_model : " << i.experiment_instrument_model << endl;
				cout << "\tsample_accession : " << accession_to_str(i.sample_accession) << endl;
				cout << "\tsample_taxa : " << i.sample_taxa << endl;

				// Due to the large amount of data involved, sample attributes are handled separately
				MAP<SraAccession, char*>::const_iterator iter = sample_attributes.find(i.sample_accession);

				if( iter != sample_attributes.end() ){

					const MAP<string, string> local = buffer_to_map(iter->second);

					for(MAP<string, string>::const_iterator j = local.begin();j != local.end();++j){
						i.sample_attributes[j->first] = j->second;
					}

				}

				if( !i.sample_attributes.empty() ){

					cout << "\tsample_attributes" << endl;

					for(MAP<string, string>::const_iterator j = i.sample_attributes.begin();
						j != i.sample_attributes.end();++j){
						
						cout << "\t\t" << j->first << " : " << j->second << endl;
					}
				}

				cout << "\tstudy_accession : " << accession_to_str(i.study_accession) << endl;
				cout << "\tstudy_title : " << i.study_title << endl;
				cout << "\tstudy_abstract : " << i.study_abstract << endl;
			}

			return EXIT_SUCCESS;
		}
		
		cerr << "Writing binary file of metadata to " << opt.output_file << " ... ";

		// Write a file of SRA metadata

		#ifdef CSV_OUTPUT
		// Write a CSV file of some of the SRA metadata
		for(deque<FilterInfo>::const_iterator i = db.begin();i != db.end();++i){
			fout << accession_to_str(i->run_accession) << ',' << i->date_received << ',' 
				<< i->experiment_library_strategy << ',' << i->experiment_library_source
				<< ',' << i->number_of_bases << endl;
		}
		#else // BINARY_OUTPUT

		// The binary output start with the number of FilterInfo elements
		binary_write( fout, db.size() );

		size_t num_sample_attribute_inject = 0;

		while( !db.empty() ){
			
			FilterInfo i = db.front();

			db.pop_front();

			MAP<SraAccession, char*>::const_iterator iter = sample_attributes.find(i.sample_accession);

			// Inject the sample attributes if present
			if( iter != sample_attributes.end() ){

				const MAP<string, string> local = buffer_to_map(iter->second);;

				for(MAP<string, string>::const_iterator j = local.begin();j != local.end();++j){
					i.sample_attributes[j->first] = j->second;
				}
				
				++num_sample_attribute_inject;
			}

			binary_write(fout, i);
		}

		#endif // CSV_OUTPUT

		// Clean up
		for(MAP<SraAccession, char*>::iterator i = sample_attributes.begin();i != sample_attributes.end();++i){

			delete [] i->second;
			i->second = NULL;
		}

		cerr << "done." << endl;
		cerr << "Injected sample attribute data for " << num_sample_attribute_inject << " SRA runs" << endl;

		profile = time(NULL) - profile;
		
		cerr << "Finished inventory in " << profile << " sec" << endl;
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

void parse_sra_text(deque<FilterInfo> &m_db, const string &m_metadata_file)
{
	// Parse the tab-delimited SRA_Accessions.tab file
	bool first_time_sra_accession_file = true;

	size_t num_sra_accession_col = 0;

	int accession_col = -1;
	int status_col = -1;
	int upate_col = -1;
	int publish_col = -1;
	int received_col = -1;
	int type_col = -1;
	int visibility_col = -1;
	int experiment_col = -1;
	int sample_col = -1;
	int study_col = -1;
	int center_col = -1;
	int spots_col = -1;
	int bases_col = -1;

	TarIterator iter(m_metadata_file);

	while(iter){

		const SRAFileType f = sra_file_type( iter.filename() );

		if(f == SRA_ACCESSIONS){

			// We expect to read a tab-delimited set of columns:
			// 0: Accession
			// 1: Submission
			// 2: Status
			// 3: Updated
			// 4: Published
			// 5: Received
			// 6: Type
			// 7: Center
			// 8: Visibility
			// 9: Alias
			// 10: Experiment
			// 11: Sample
			// 12: Study
			// 13: Loaded
			// 14: Spots
			// 15: Bases
			// 16: Md5sum
			// 17: BioSample
			// 18: BioProject
			// 19: ReplacedBy

			// Dates are formatted as YYYY-MM-DDThh:mm:ssZ
			// 	For example: "2010-03-24T03:10:22Z"

			if(first_time_sra_accession_file){

				const deque<string> cols = split(*iter, '\t');

				num_sra_accession_col = cols.size();

				first_time_sra_accession_file = false;

				// Reset the column values to invalid values
				// (since there are multiple tab-delimited files
				// that share some of the column headers)

				accession_col = -1;
				status_col = -1;
				upate_col = -1;
				publish_col = -1;
				received_col = -1;
				type_col = -1;
				visibility_col = -1;
				experiment_col = -1;
				sample_col = -1;
				study_col = -1;
				center_col = -1;
				spots_col = -1;
				bases_col = -1;

				for(deque<string>::const_iterator i = cols.begin();i != cols.end();++i){

					if(*i == "Accession"){
						accession_col = i - cols.begin();
					}

					if(*i == "Status"){
						status_col = i - cols.begin();
					}

					if(*i == "Updated"){
						upate_col = i - cols.begin();
					}

					if(*i == "Published"){
						publish_col = i - cols.begin();
					}

					if(*i == "Received"){
						received_col = i - cols.begin();
					}

					if(*i == "Type"){
						type_col = i - cols.begin();
					}

					if(*i == "Visibility"){
						visibility_col = i - cols.begin();
					}

					if(*i == "Experiment"){
						experiment_col = i - cols.begin();
					}

					if(*i == "Sample"){
						sample_col = i - cols.begin();
					}

					if(*i == "Study"){
						study_col = i - cols.begin();
					}

					if(*i == "Center"){
						center_col = i - cols.begin();
					}

					if(*i == "Spots"){
						spots_col = i - cols.begin();
					}

					if(*i == "Bases"){
						bases_col = i - cols.begin();
					}
				}

				// Make sure that we found all of the expected columns
				if(accession_col < 0){
					throw __FILE__ ":parse_sra_text: Did not find \"Accession\" column in SRA Accessions file";
				}

				if(status_col < 0){
					throw __FILE__ ":parse_sra_text: Did not find \"Status\" column in SRA Accessions file";
				}

				if(upate_col < 0){
					throw __FILE__ ":parse_sra_text: Did not find \"Updated\" column in SRA Accessions file";
				}

				if(publish_col < 0){
					throw __FILE__ ":parse_sra_text: Did not find \"Published\" column in SRA Accessions file";
				}

				if(received_col < 0){
					throw __FILE__ ":parse_sra_text: Did not find \"Received\" column in SRA Accessions file";
				}

				if(type_col < 0){
					throw __FILE__ ":parse_sra_text: Did not find \"Type\" column in SRA Accessions file";
				}

				if(visibility_col < 0){
					throw __FILE__ ":parse_sra_text: Did not find \"Visibility\" column in SRA Accessions file";
				}

				if(experiment_col < 0){
					throw __FILE__ ":parse_sra_text: Did not find \"Experiment\" column in SRA Accessions file";
				}

				if(sample_col < 0){
					throw __FILE__ ":parse_sra_text: Did not find \"Sample\" column in SRA Accessions file";
				}

				if(study_col < 0){
					throw __FILE__ ":parse_sra_text: Did not find \"Study\" column in SRA Accessions file";
				}

				if(center_col < 0){
					throw __FILE__ ":parse_sra_text: Did not find \"Center\" column in SRA Accessions file";
				}

				if(spots_col < 0){
					throw __FILE__ ":parse_sra_text: Did not find \"Spots\" column in SRA Accessions file";
				}

				if(bases_col < 0){
					throw __FILE__ ":parse_sra_text: Did not find \"Bases\" column in SRA Accessions file";
				}
			}
			else{

				const deque<string> cols = split(*iter, '\t');

				if(cols.size() != num_sra_accession_col){
					throw __FILE__ ":parse_sra_text: Did not read the expected number of columns in the SRA Accession file";
				}

				// Only check RUN accessions
				if(cols[type_col] == "RUN"){

					const SraAccession local_accession = str_to_accession(cols[accession_col]);

					bool valid = true;

					if( ( cols[status_col] == "suppressed" ) || 
						( cols[status_col] == "controlled_access" ) || 
						( cols[status_col] == "unpublished" ) || 
						( cols[visibility_col] == "suppressed" ) ||
						( cols[visibility_col] == "controlled_access" ) ){
						
						valid = false;
					} 

					// Only store valid SRA runs
					if(valid){

						// Lookup or create an entry in our database
						m_db.push_back( FilterInfo() );

						FilterInfo &ref = m_db.back();

						ref.run_accession = local_accession;
						
						ref.valid = valid;

						if(cols[spots_col] != "-"){
							ref.number_of_spots = str_to_size_t(cols[spots_col]);
						}

						if(cols[bases_col] != "-"){
							ref.number_of_bases = str_to_size_t(cols[bases_col]);
						}

						ref.date_received = Date(cols[received_col]);

						if(cols[experiment_col] != "-"){
							ref.experiment_accession = str_to_accession(cols[experiment_col]);
						}

						if( (cols[sample_col] != "-") && (cols[sample_col] != "Multiplex") ){
							ref.sample_accession = str_to_accession(cols[sample_col]);
						}

						if(cols[study_col] != "-"){
							ref.study_accession = str_to_accession(cols[study_col]);
						}

						if(cols[center_col] != "-"){
							ref.sample_attributes["Center"]  = cols[center_col];
						}
					}
				}
			}
		}

		++iter;
	}
}

void parse_sra_metadata(deque<FilterInfo> &m_db, MAP<SraAccession, char*> &m_sample_attributes,
	const string &m_metadata_file,
	const InventoryOptions &m_opt)
{
	cerr << "Parsing the tab-delimited tables ... ";

	// Step 1: Parse the tab-delimited text sub-files in the metadata file
	parse_sra_text(m_db, m_metadata_file);

	cerr << "found " << m_db.size() << " SRA runs" << endl;

	cerr << "Parsing the XML data ... ";

	// Step 2: Parse the xml files
	MAP<SraAccession, FilterInfo> xml_info;

	// Periodically purge the number XML records after updating the 
	// final database
	const size_t max_num_xml = 100000;

	string curr_filename = "";

	SraAccession experiment_accession = INVALID_ACCESSION;
	SraAccession sample_accession = INVALID_ACCESSION;
	SraAccession study_accession = INVALID_ACCESSION;
	bool sample_attribute = false;
	string sample_attribute_tag;
	string sample_attribute_value;

	size_t num_experiment = 0;
	size_t num_sample = 0;
	size_t num_study = 0;

	TarIterator iter(m_metadata_file);

	while(iter){

		if( curr_filename != iter.filename() ){

			// Reset the state variables
			experiment_accession = INVALID_ACCESSION;
			sample_accession = INVALID_ACCESSION;
			study_accession = INVALID_ACCESSION;
			sample_attribute = false;
			sample_attribute_tag.clear();
			sample_attribute_value.clear();

			curr_filename = iter.filename();

			// Limit the number of XML records to store RAM (to avoid running out of
			// memory!)
			if(xml_info.size() >= max_num_xml){

				merge_db(m_db, m_sample_attributes, xml_info, num_experiment, num_sample, num_study);

				xml_info.clear();
			}
		}

		const SRAFileType f = sra_file_type( iter.filename() );

		if(f == EXPERIMENT_XML){
				
			// Append a space to the search key "EXPERIMENT " to avoid picking up
			// "EXPERIMENT_SET".
			if( (*iter).find("<EXPERIMENT ") != string::npos ){
				experiment_accession = str_to_accession( parse_key(*iter, "accession") );
			}
				
			if( (*iter).find("<TITLE>") != string::npos){
				
				if( experiment_accession == INVALID_ACCESSION){
					throw __FILE__ ":parse_sra_metadata: Orphaned experiment title";
				}
				
				xml_info[experiment_accession].experiment_title = 
					parse_xml("TITLE", *iter);
			}
				
			if( (*iter).find("<DESIGN_DESCRIPTION>") != string::npos){
				
				if(experiment_accession == INVALID_ACCESSION){
					throw __FILE__ ":parse_sra_metadata: Orphaned experiment design description";
				}
				
				xml_info[experiment_accession].experiment_design_description = 
					parse_xml("DESIGN_DESCRIPTION", *iter);
			}

			if( (*iter).find("<LIBRARY_NAME>") != string::npos){
				
				if(experiment_accession == INVALID_ACCESSION){
					throw __FILE__ ":parse_sra_metadata: Orphaned experiment library name";
				}
				
				xml_info[experiment_accession].experiment_library_name = 
					parse_xml("LIBRARY_NAME", *iter);
			}

			if( (*iter).find("<LIBRARY_STRATEGY>") != string::npos){
				
				if(experiment_accession == INVALID_ACCESSION){
					throw __FILE__ ":parse_sra_metadata: Orphaned experiment library strategy";
				}
				
				xml_info[experiment_accession].experiment_library_strategy = 
					parse_xml("LIBRARY_STRATEGY", *iter);
			}

			if( (*iter).find("<LIBRARY_SOURCE>") != string::npos){
				
				if(experiment_accession == INVALID_ACCESSION){
					throw __FILE__ ":parse_sra_metadata: Orphaned experiment library source";
				}
				
				xml_info[experiment_accession].experiment_library_source = 
					parse_xml("LIBRARY_SOURCE", *iter);
			}

			if( (*iter).find("<LIBRARY_SELECTION>") != string::npos){
				
				if(experiment_accession == INVALID_ACCESSION){
					throw __FILE__ ":parse_sra_metadata: Orphaned experiment library selection";
				}
				
				xml_info[experiment_accession].experiment_library_selection = 
					parse_xml("LIBRARY_SELECTION", *iter);
			}

			if( (*iter).find("<INSTRUMENT_MODEL>") != string::npos){
				
				if(experiment_accession == INVALID_ACCESSION){
					throw __FILE__ ":parse_sra_metadata: Orphaned experiment instrument model";
				}
				
				xml_info[experiment_accession].experiment_instrument_model = 
					parse_xml("INSTRUMENT_MODEL", *iter);
			}

			// Check for the presence of dbgap data for which access is controlled.
			// See: https://www.ncbi.nlm.nih.gov/books/NBK56558/
			if( (*iter).find("<EXTERNAL_ID namespace=\"dbgap\">") != string::npos){
				
				if(experiment_accession == INVALID_ACCESSION){
					throw __FILE__ ":parse_sra_metadata: Orphaned experiment dbgap id";
				}
				
				xml_info[experiment_accession].valid = false;
			}
		}

		if(f == SAMPLE_XML){

			if( (*iter).find("<SAMPLE ") != string::npos ){
				sample_accession = str_to_accession( parse_key(*iter, "accession") );
			}
			
			if( (*iter).find("<SCIENTIFIC_NAME>") != string::npos){
				
				if(sample_accession == INVALID_ACCESSION){
					throw __FILE__ ":parse_sra_metadata: Orphaned sample scientific name";
				}
				
				xml_info[sample_accession].sample_taxa = 
					parse_xml("SCIENTIFIC_NAME", *iter);
			}
				
			if( (*iter).find("<SAMPLE_ATTRIBUTE>") != string::npos){
				sample_attribute = true;
			}
				
			if( (*iter).find("</SAMPLE_ATTRIBUTE>") != string::npos){
				sample_attribute = false;
			}
				
			if( sample_attribute && ( (*iter).find("<TAG>") != string::npos ) ){
				sample_attribute_tag = parse_xml("TAG", *iter);
			}
				
			if( sample_attribute && ( (*iter).find("<VALUE>") != string::npos ) ){
			
				sample_attribute_value = parse_xml("VALUE", *iter);
				
				if( sample_attribute_tag.empty() && (sample_accession == INVALID_ACCESSION) ){
					throw __FILE__ ":parse_sra_metadata: Orphaned sample attribute value";
				}
				
				// Skip the BioSampleModel tagged values (these entries are *not* reported
				// on the SRA BioSample web-pages).
				if(sample_attribute_tag != "BioSampleModel"){
				
					xml_info[sample_accession].sample_attributes[sample_attribute_tag] = 
						sample_attribute_value;
				}
			}
		}
		
		if(f == STUDY_XML){

			// Add a space after "<STUDY" to avoid picking up "<STUDY_SET"
			if( (*iter).find("<STUDY ") != string::npos){
				study_accession = str_to_accession( parse_key(*iter, "accession") );
			}

			if( (*iter).find("<STUDY_TITLE>") != string::npos){
				
				if(study_accession == INVALID_ACCESSION){

					cerr << "XML: " << iter.filename() << endl;
					throw __FILE__ ":parse_sra_metadata: Orphaned study title";
				}

				xml_info[study_accession].study_title = parse_xml("STUDY_TITLE", *iter);
			}

			if( (*iter).find("<STUDY_ABSTRACT>") != string::npos ){
				
				if(study_accession == INVALID_ACCESSION){
					throw __FILE__ ":parse_sra_metadata: Orphaned study abstract";
				}

				xml_info[study_accession].study_abstract = parse_xml("STUDY_ABSTRACT", *iter);
			}
		}

		++iter;
	}

	// Perform a final merge of the xml-derived information
	merge_db(m_db, m_sample_attributes, xml_info, num_experiment, num_sample, num_study);

	cerr << "done." << endl;

	const size_t num_sra = m_db.size();

	if(num_sra > 0){

		cerr << "Found XML annotation for:" << endl;
		cerr << '\t'<< num_experiment << " (" << (100.0*num_experiment)/num_sra 
			<< "%) SRA runs by association with SRA experiments" << endl;
		cerr << '\t'<< num_sample << " (" << (100.0*num_sample)/num_sra 
			<< "%) SRA runs by association with SRA samples" << endl;
			cerr << '\t'<< m_sample_attributes.size() << " (" 
				<< (100.0*m_sample_attributes.size())/num_sra 
			<< "%) SRA sample records have attribute data (to be added later)" << endl;
		cerr << '\t'<< num_study << " (" << (100.0*num_study)/num_sra 
			<< "%) SRA runs by association with SRA studies" << endl;
	}
}


void merge_db(deque<FilterInfo> &m_db, MAP<SraAccession, char*> &m_sample_attributes,
	const MAP<SraAccession, FilterInfo> &m_xml_info, 
	size_t &m_num_experiment, size_t &m_num_sample, size_t &m_num_study)
{
	for(deque<FilterInfo>::iterator r = m_db.begin();r != m_db.end();++r){

		if(r->experiment_accession != INVALID_ACCESSION){

			MAP<SraAccession, FilterInfo>::const_iterator x = 
				m_xml_info.find(r->experiment_accession);

			if( x != m_xml_info.end() ){
				
				bool updated = false;

				if( !x->second.valid && r->valid){

					r->valid = false;
					updated = true;
				}

				#define	UPDATE(VAR) \
					if( !x->second.VAR.empty() ){\
						r->VAR = x->second.VAR;\
						updated = true;\
					}

				UPDATE(experiment_title)
				UPDATE(experiment_library_name)
				UPDATE(experiment_library_strategy)
				UPDATE(experiment_library_source)
				UPDATE(experiment_library_selection)
				UPDATE(experiment_instrument_model)

				#undef UPDATE

				if(updated){
					++m_num_experiment;
				}
			}
		}

		if(r->sample_accession != INVALID_ACCESSION){

			MAP<SraAccession, FilterInfo>::const_iterator s = 
				m_xml_info.find(r->sample_accession);

			if( s != m_xml_info.end() ){

				bool updated = false;

				#define	UPDATE(VAR) \
					if( !s->second.VAR.empty() ){\
						r->VAR = s->second.VAR;\
						updated = true;\
					}

				UPDATE(sample_taxa)

				// Are there attributes to store?
				if( !s->second.sample_attributes.empty() ){

					MAP<SraAccession, char*>::iterator attrib_iter = 
						m_sample_attributes.find(r->sample_accession);
					
					MAP<string, string> local;

					if( ( attrib_iter != m_sample_attributes.end() ) && 
						(attrib_iter->second != NULL) ){

						local = buffer_to_map(attrib_iter->second);

						// Remove the existing buffer
						delete [] attrib_iter->second;
						attrib_iter->second = NULL;
					}

					// Add the new sample attributes to any existing attributes
					for(MAP<string, string>::const_iterator i = s->second.sample_attributes.begin();
						i != s->second.sample_attributes.end();++i){

							local[i->first] = i->second;
					}

					// Store the new buffer
					m_sample_attributes[r->sample_accession] = map_to_buffer(local);
				}

				//if( !s->second.sample_attributes.empty() ){
				//	
				//	for(MAP<string, string>::const_iterator i = s->second.sample_attributes.begin();
				//		i != s->second.sample_attributes.end();++i){
				//
				//			r->sample_attributes[i->first] = i->second;
				//			updated = true;
				//	}
				//}

				#undef UPDATE

				if(updated){
					++m_num_sample;
				}
			}
		}

		if(r->study_accession != INVALID_ACCESSION){

			MAP<SraAccession, FilterInfo>::const_iterator s = 
				m_xml_info.find(r->study_accession);

			if( s != m_xml_info.end() ){

				bool updated = false;

				#define	UPDATE(VAR) \
					if( !s->second.VAR.empty() ){\
						r->VAR = s->second.VAR;\
						updated = true;\
					}
				
				UPDATE(study_title)
				UPDATE(study_abstract)

				if(updated){
					++m_num_study;
				}
			}
		}
	}
}

inline bool find_extension(const string &m_filename, const string &m_ext)
{
	const string::size_type loc = m_filename.rfind(m_ext);
	
	if(loc == string::npos){
		return false;
	}
	
	return ( ( loc + m_ext.size() ) == m_filename.size() );
}

// Match the file type by looking at the file name suffix
SRAFileType sra_file_type(const string &m_filename)
{
	
	if( find_extension(m_filename, ".run.xml") ){
		return RUN_XML;
	}
	
	if( find_extension(m_filename, ".experiment.xml") ){
		return EXPERIMENT_XML;
	}
	
	if( find_extension(m_filename, ".sample.xml") ){
		return SAMPLE_XML;
	}
	
	if( find_extension(m_filename, ".study.xml") ){
		return STUDY_XML;
	}
	
	if( find_extension(m_filename, "SRA_Accessions") ){
		return SRA_ACCESSIONS;
	}
	
	if( find_extension(m_filename, "SRA_Run_Members") ){
		return SRA_RUN_MEMBERS;
	}

	return UNKNOWN;
}

string parse_xml(const string &m_key, const string &m_buffer)
{
	string query = "<" + m_key + ">";
	
	size_t begin = m_buffer.find(query);
	
	if(begin == string::npos){
		throw __FILE__ ":parse_xml: Unable to find begin key";
	}
	
	begin += query.size();
	
	query = "</" + m_key + ">";
	
	const size_t end = m_buffer.rfind(query);
	
	if(end == string::npos){
		throw __FILE__ ":parse_xml: Unable to find end key";
	}
	
	if(begin > end){
		throw __FILE__ ":parse_xml: XML delimeters are out of order";
	}
	
	return m_buffer.substr(begin, end - begin);
}

string parse_key(const string &m_buffer, string m_key /* make a copy */)
{
	const size_t len = m_buffer.size();
	
	// Add the key delimiter
	m_key += "=";

	size_t loc = m_buffer.find(m_key);
	
	if(loc == string::npos){

		//cerr << "m_buffer = " << m_buffer << endl;
		//cerr << "m_key = " << m_key << endl;
		throw __FILE__ ":parse_key: Unable to find key=";
	}
	
	loc += m_key.size() + 1; // Add one to skip the opening " symbol
	
	if(loc >= len){
		throw __FILE__ ":parse_key: Missing key value";
	}

	for(size_t i = 1;(loc + i) < len;++i){
		
		if(m_buffer[loc + i] == '"'){
			return m_buffer.substr(loc, i);
		}
	}
	
	throw __FILE__ ":parse_key: Unable to extract key -- no closing quote";
	
	return string();
}

MAP<string, string> buffer_to_map(char *m_buffer)
{
	MAP<string, string> ret;

	char* key_begin = m_buffer;

	while(true){

		// Stop reading the buffer when we find an end-of-string symbol (\0) in
		// the key position.
		if(*key_begin == '\0'){
			break;
		}

		char* key_end = key_begin;

		while(*key_end != '\0'){
			++key_end;
		}

		char* value_begin = key_end + 1;

		char* value_end = value_begin;

		while(*value_end != '\0'){
			++value_end;
		}

		ret[string(key_begin, key_end)] = string(value_begin, value_end);

		key_begin = value_end + 1;
	}

	return ret;
}

char* map_to_buffer(const MAP<string, string> &m_map)
{
	size_t len = 1; // Must store an extra '\0' to terminate the buffer

	for(MAP<string, string>::const_iterator i = m_map.begin();i != m_map.end();++i){
		len += i->first.size() + i->second.size() + 2; // Store the end of string delimeters
	}

	char *buffer = new char[len];

	if(buffer == NULL){
		throw __FILE__ ":map_to_buffer: Unable to allocate buffer";
	}

	char *ptr = buffer;

	for(MAP<string, string>::const_iterator i = m_map.begin();i != m_map.end();++i){

		memcpy(ptr, i->first.c_str(), i->first.size() + 1);
		ptr += i->first.size() + 1;

		memcpy(ptr, i->second.c_str(), i->second.size() + 1);
		ptr += i->second.size() + 1;
	}

	// Terminate the buffer with a '\0' in the key (i.e. first) position
	*ptr = '\0';

	return buffer;
}
