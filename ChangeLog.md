# Maestro changes
- **Version 0.9** October 19, 2020
	- Improved the buffering of the transposed bitmatrix in build_db() to significantly reduce the number
	  write operations (and hopefully fix the issues with Lustre/FSX stalling/hanging ...).

- **Version 0.8** September 13, 2020
	- Provide for separate scratch paths: `--scratch.bloom` and `--scratch.database`
	- Fixed template specialization-related compiler errors on g++
	- Added the option to `--skip` user-specified SRA run accessions. This is useful when trying to
	  deal with accessions that cannot be downloaded and/or parsed by the SRA toolkit (i.e. ERR1197571).
	  	- At some point, we will need to add an "include" option to be able process previously skipped
		  accessions (should they get fixed and become possible to include again).
	- Fixed bug in final database creation step (when all SRA records have been processed)

- **Version 0.7** August 29, 2020
	- Track worker memory usage to check for memory leaks
	- Even when streaming SRA records, remove *.cache and other SRA record-specific files.

- **Version 0.6** August 15, 2020
	- Handle aligned colorspace SRA records as a special case. These records cannot be read using the
	strategy of loading primary alignments followed by unaligned reads.
	- Use the SRA metadata to extract the number of bases and scale the size of the counting Bloom
	filters (to obtain the desired false positive rate).

- **Version 0.5** August 5, 2020
	- Added code to track the:
		- "deflation"; the number of Bloom filter bytes/number of sequence bytes
		- rate of SRA processing in kmers/second and bp/sec
		- ratio of kmers/bp for each Bloom filter -- what fraction of sequence is unique?
	- Modified the process_event() loop to automatically retry failed download and Bloom filter attempts
	  (up to the number of allowed retries).
	- Added an option to force a retry of (STATUS_BLOOM_FAILURE, i.e. "hard") failed Bloom filters. 
	  Use "--retry.bloom" to retry all failed Bloom filters.
	- Added a new command line option, "--delay <number of seconds>", to ensure that download and/or streaming requests
	  do not occur more frequently than once per the number of specified seconds.
	- Added additional reporting to the make_bloom_filter() function to track the progress reading through the SRA
	  record. This will hopefully inform on the utility of restarting streaming failures during Bloom filter construction.
	- Catch SRA NGS error messages to help diagnose connectivity issues

- **Version 0.4** July 30, 2020
	- Modified the restore_bloom() in maestro_main.cpp to restore SRA accessions that are labeled
	  as STATUS_DATABASE_FAIL (in addition to STATUS_BLOOM_SUCCESS). This will recover from previous database creation failures.
	- Moved some file I/O code to a new source file (file_io.cpp) from the maestro_main.cpp file.
	- Created a new help program, "manual_db", to update the accession in a status file with the accessions
	  in a database file that is being manually copied to S3. This is needed when we get an upload failure
	  (due to "aws s3 mv/cp" failing to upload a database file -- most likely due to Lustre FSX being overtaxed
	  during the creation of a large databaes file.)

# SRA inventory changes
- **Version 0.7** September 24, 2020
	- Added the option to specify an optional list of SRA run accessions to include.

# CALDERA changes
- **Version 0.4c** December 10, 2020
	- Changed program name from `bigsi++` to `caldera`.
- **Version 0.4b** September 25, 2020
	- Fixed supurious output when the input query is too small for a single kmer.
