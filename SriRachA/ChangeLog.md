# SriRachA
A parallel tool for k-mer based searchig of SRA records

## Change Log
	-Version 0.43 (June 1, 2020)
		- Added the ability to search local sra files.
		- Improved match threshold test to remove floating point truncation issues.
	-Version 0.42 (May 14, 2020)
		- Added the sequence of the matching reads to the output file.
	-Version 0.41 (May 13, 2020)
		- Added a symbol to the output (i.e. "//") to indicate that all SRA accessions were successfully searched.
	-Version 0.4 (May 10, 2020)
		- Added additional retry attempts (with sleeping) to the sra_stream code to handle cases when the 
		  VDB API appears to give up too easily (in calls to VCursorCellDataDirect and VDBManagerPathType).
		- Fixed the base counting in the search statistics.
		- Compared two strategies for computing k-mer set overlap: (a) lookup the read k-mers in the set of
		  target k-mers using lower_bound(), and (b) iterating through both sorted lists of k-mers and
		  counting the matching k-mers. Since the read k-mer set is much small, the lower_bound() strategy
		  is significantly faster (it does not clobber the cache).
	-Version 0.3 (May 4, 2020)
		- (Hopefully) fixed the problem in sra_stream.cpp that, on rare occasions, the string length returned by 
		  VCursorCellDataDirect does not equal the sum of read lengths (also returned by VCursorCellDataDirect).
		- Added reporting of SRA search statistics: number of reads and number of bases searched.
	
