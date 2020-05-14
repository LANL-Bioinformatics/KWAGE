# SriRachA: Per-read, streaming searches of the SRA (Version 0.42, May 14, 2020)

The SriRachA program queries all of the reads in a set of user specified SRA accessions against a user-provided set of query sequences. This capability is useful for confirming and investigating the per-accession matches identified by the BIGSI++ search tool.

- Searches are kmer-based:
  - Each query sequence is converted into a set of canonical kmers
  - Each read in every SRA accession is converted into a set of canonical kmers.
  - Matching reads are identified when the fraction of read kmers that are found in the query kmer set exceeds a user-specified threshold.
- The SRA data is streamed using the VDB API provided as part of the SRA tools
  - To accelerate searches, the Message Passing Interface (MPI) is used to search non-overlapping slices of the sequence read data in an SRA accession. This can increase the bandwidth when downloading data from the SRA.
  - Currenlty, only the sequence data from an SRA record is downloaded (i.e. no quality scores, read names or other metadata is downloaded).

Due to the large amounts of data that may need to be downloaded, it is *strongly* recommended to run SriRachA in the cloud, using either Amazon (AWS) or Google (GCP). NCBI has placed in the SRA data in both the AWS and Google clouds, and you should make sure to run in the same region as the data (currently `us-east-1` for AWS and any region that starts with `us-` for GCP).

Attempting to stream SRA data to a private (non-cloud) network will likely be painfully slow!

## Building SriRachA from source

The SriRachA program is written in C++ and depends the SRA tools provided by the NCBI. Please see the primary BIGSI++ [documentation](https://github.com/LANL-Bioinformatics/BIGSI-plus-plus/blob/master/README.md) for instructions on how to build and install the SRA tools. MPI is used to search SRA records in parallel, so you will need to have an MPI development kit installed.

There is a simple `Makefile` that will need to be manually edited to compile SriRachA. Please set the `SRA` variable to the location where the SRA tools are installed. By default, this variable is equal to `$(HOME)/src/BIGSI/SRA`, but can be any directory. OpenMP is *not* currently used, but this may change in the near future...

Type `make` to compile the program. This should create the executable `sriracha`. There is no `make install`, but you can copy this executable file to the location of your choice.

If you run into any problems compiling, please let me know.

## Running SriRachA

Running `sriracha` without any command line arguments will display the list of allowed arguments:

```
Usage for SriRachA (v. 0.4):
	-i <input sequence files> (can be repeated)
	[-o <output filename>] (default is stdout)
	[--read.len.min <minimum read length>] (default is 0)
	[--max-results <maximum number of results to show per accession/query>] (default is 100)
	[-a <list of SRA accessions in a text file>]
	[-v (increase the verbosity: silent, tacitern, normal, chatty. Default is silent)]
	[--retry <maximum number of download atttemps>] (default is 0)
	Search strategies
		[--search-by-align]
			(Not implemented yet!)
		[--search-by-kmer] (default)
			[-k <k-mer length>] (default is 11)
			[-t <match threshold>] (default is 0.8)
			[-n <min number valid kmer>] (default is 1)
			[--read.complexity.min <min read complexity>] (default is 0.75)
		[--search-by-bloom]
			(Not implemented yet!)
	<SRA accession1> ...
```

- The `-i` flag specifies the input file containing a (set) of nucleotide sequences in either FASTA or FASTQ formats. Each sequence is treated as a separate query. Multiple queries are allowed (either within a file, or by providing multiple files).
- The optional `-o` flag specifies the output file. If no output file is provided, the output is written to stdout.
- The optional `--read.len.min` flag specifies the minimum length of an SRA read. Reads with fewer bases will not be searched.
- The optional `--max-results` flag specifies the maximum number of matching reads to report for each SRA accession/query pair. This is intended to keep the output files from being swamped in the cases where there are millions of reads that match the specified query sequences.
- The optional `-a` flag allows the set of SRA accessions to be searched to be specified in a file (one accession per line), in addition to any accession that are specified via the command line.
- The optional `-v` flag controls the amount of information that is provided to the user while the program is running. Be warned that "chatty" is intended for debugging, and will produce copious amounts of output!
- Finally, specify any number of SRA accessions to search on the command line (without any flags). If no SRA accessions are specified on the command line or in a file (via the `-a` flag), then the program will read SRA accession from stdin (for integration with other tools via pipes).

Only one search stragtegy is currently implemented: `search-by-kmer`
- The optional `-k` flag controls the kmer length. The kmer length must be greater than, or equal to, 3 and less than, or equal to, 32.
- The optional `-t` flag specifies the fraction of read kmers that must match a query kmer to report a match
- The optional `-n` flag specifies the minimum number of kmers that each read must contain to be searched. This is intended to allow short reads to be skipped.
- The optional `--read.complexity.min` flag specifies the minimim complexity that a read must have before it is searched. The read complexity is defined as: (number of unique kmers in each read)/(total number of kmers in each read).

### Output format

Only a single output format is currently provided. The output of SriRachA is five, ASCII-text, tab-separated columns:
- The first column is the SRA accession
- The second column is the `<read index>.<read pair id>`. The read index is the same read index used by the SRA database, to facilitate the extraction of the actual read sequence/quality from the [SRA database](https://trace.ncbi.nlm.nih.gov/Traces/sra/).
- The third column is the match score, defined as (number of shared read/query kmers)/(number of unique read kmers).
- The fourth column is the read sequence.
- The fifth column is the name of the matching query sequence (extracted from the FASTA/FASTQ defline).

Finally, the end of the file is signified by the string `//`. This indicates when the search has been completed, and is useful when many thousands of SRA accessions are being queried in parallel.
