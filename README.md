# CALDERA: Compressed Approach for Low-overhead Digital Exploration of Read Archives

A C++ implementation of the ultra-fast, [Bloom filter](https://en.wikipedia.org/wiki/Bloom_filter)-based sequence search originally developed by the [Iqbal group](https://www.nature.com/articles/s41587-018-0010-1). This project is similar in spirit to the Iqbal group's BIGSI follow-on project: [COBS: a Compact Bit-Sliced Signature Index](https://arxiv.org/abs/1905.09624). While the goals are very similar, the CALDERA project is an independent implementation of the original BIGSI algorithm with the addition of the following features:
1. Parallel, C++ implementation (using both OpenMP and MPI)
2. Adaptive Bloom filter size (similar to COBS)
3. Compressed database files
4. Bloom filter construction directly from [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) (SRA) files

The CALDERA pipeline has several components:
1. `sra_inventory`
	- Parse and convert the NCBI-supplied [SRA metadata information](tp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata) into a single binary file. This file is read by the `maestro` program to add metadata and annotation information to the final database files.
2. `maestro`
	- Streaming Bloom filter construction from SRA accessions.
	- Aggregation and transposition of of individual Bloom filters into database files.
3. `caldera`
	- Searching the resulting database with a nucleic acid query:.

The following helper applications are also provided:
1. `dump_db`, a tool for dumping the header, annotation, and first few bit slices of a database file.
2. `SriRachA`, a tool for per-read, kmer-based searching of SRA records against a set of user-provided query sequences. This tool is useful for confirming and investigating the per-accession matches reported by CALDERA. This tool is contained in the `SriRachA/` subdirectory (with a separate `Makefile` and `documentation`).

## Building and installing the code
The CALDERA code is written in C++ and requires the MPI (message passing interface; the code has been tested using [OpenMPI](https://www.open-mpi.org/)). A compiler that supports OpenMP is also suggested, but not required.

Currently, downloading the SRA records is orchestrated by the `maestro` program. While previous versions relied on a cluster scheduling system (like Slurm), this dependancy has been **removed**. As discussed below, downloading the entire SRA to a private network (i.e. non-[AWS](https://aws.amazon.com/), non-[Google](https://cloud.google.com/)) is **not** going to work -- downloading will take too long for the several petabytes of SRA data.

During the database creation process, database files are written to an AWS S3 bucket (to avoid exhausting expensive local storage). AWS command line tools are used to copy database files to S3. **This storage approach will need to be modified to run on a different cloud provider (i.e. Google).** The dependancy on AWS command line tools can be avoided by specifying the `--s3.no-write` flag to the `maestro` program. Depending on the number of SRA records being indexed, the storage requirements for the database files can be significant!

### Step 1: Download the SRA toolkit and associated APIs
In addition, the SRA toolkit and associated APIs are required if you would like to be able to download and/or construct Bloom filters from SRA files. In particular, you will need to download and build the following three GitHub projects from NCBI:
1. [ncbi-vdb](https://github.com/ncbi/ncbi-vdb)
2. [ngs](https://github.com/ncbi/ngs)
3. [sra-tools](https://github.com/ncbi/sra-tools)

Please note that these NCBI libraries are only needed if you will need to either download SRA files and/or create Bloom filters from SRA files. If you only need to either build a database from existing Bloom filters or search an existing database, then these NCBI libraries are *not* needed.

#### How to actually build the SRA toolkit packages
It can be a little tricky to compile the ncbi-vdb, ngs and sra-tools packages. Here is a step-by-step guide (with specific package install commands for a Red Hat-like OS that uses the `yum` package manager):

1. Install the `git` repository control software (e.g. `sudo yum install git`).
2. Install a c/c++ compiler (e.g. `sudo yum group install "Development Tools"`).
4. Install a Java development kit (while Java is not used in CALDERA, a jdk is needed to compile the ncbi.ngs package).
5. Install libxml2 development libraries (e.g. `sudo yum install libxml2-devel.x86_64`). Many Linux distributions ship with libxml2 but not the development headers, which will be needed.
6. Create a directory to hold all of the NCBI SRA software (e.g. `mkdir SRA`). The rest of these instructions assume that the path to the SRA packages is `$HOME/SRA`.
7. Download sra-tools from GitHub (e.g. `git clone https://github.com/ncbi/sra-tools.git`).
8. Download ngs from GitHub (e.g. `git clone https://github.com/ncbi/ngs.git`).
9. Download ncbi-vdb from GitHub (e.g. `git clone https://github.com/ncbi/ncbi-vdb.git`).
10. Compile the NCBI packages in the **following order** (order is important!):

```
# Build and install ncbi-vdb
cd $HOME/SRA/ncbi-vdb
./configure --prefix=$HOME/SRA \
	--with-ngs-sdk-prefix=$HOME/SRA/ngs

# Before building on Amazon's AWS (and depending on choice of image/instance), you may need to remove the `-Wa,-march=generic64+sse4` flags from `ncbi-vdb/libs/krypto/Makefile` prior to running `make`
make
make install
```

```
# Build and install ngs/ngs-sdk (a sub-package of ngs)
cd $HOME/SRA/ngs/ngs-sdk
./configure --prefix=$HOME/SRA/ngs/ngs-sdk
make
make install
```

```
# Build and install ngs/ngs-java
cd $HOME/SRA/ngs/ngs-java
./configure --prefix=$HOME/SRA
make
make install
```

```
# Build and install ngs
cd $HOME/SRA/ngs
./configure --prefix=$HOME/SRA \
	--with-ncbi-vdb-prefix=$HOME/SRA/ncbi-vdb \
	--with-ngs-sdk-prefix=$HOME/SRA/ngs/ngs-sdk
make
make install
```

```
# Build and install sra-tools
cd cd $HOME/SRA/sra-tools
./configure --prefix=$HOME/SRA \
	--with-ngs-sdk-prefix=$HOME/SRA/ngs/ngs-sdk \
	--with-ncbi-vdb-sources=$HOME/SRA/ncbi-vdb

make
make install
```
11. Run the vdb configuration tool: `$HOME/SRA/bin/vdb-config --interactive`.
	- Select the `[X] report cloud instance identity` option
	- At this time, do *not* select the `[ ] accept charges for AWS`.
	- Select `Save`
	- Select `Exit`

The options you select here will depend on your local compute environment (e.g. cloud vs local).

### Step 2: Edit the CALDERA Makefile
After the three separate NCBI software packages have been downloaded and compiled, you will need to edit the CALDERA Makefile to specify the locations of both the include and library directories for the ncbi-vdb and ngs packages.

Edit the `SRA_LIB_PATH` variable to specify the directory on your system that contains the SRA library files.

Edit the `SRA_INCLUDE_PATH` variable specify the directory on your system that contains the SRA include files.

Edit the `NCBI_VDB_INCLUDE_PATH` variable specify the directory on your system that contains the NCBI-VDB include files.

### Step 3: Run `make -j`

Running `make` in the BIGS++ directory will build all pipeline components. There is no separate installation step required. If you encounter any problems, please email `jgans@lanl.gov` for assistance.

## Format the NCBI SRA metadata information

The `sra_inventory` program is responsible for formatting the NCBI SRA metadata for inclusion in the final database files. The NCBI provides a compressed tar file of XML metadata dowloaded on their [FTP site](tp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata). **Do not uncompress or untar this file!** Since the uncompressed/untarred files are large and cumbersome to manage, the `sra_inventory` program *directly* reads the SRA metadata (as a compressed tar file) and converts this data into a single binary file that can be read by the `maestro` program that builds the database files from SRA records. The command line arguments for the `sra_inventory` program are:
```
Usage for SRA inventory:
	-i <XML metadata file from ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata>
	[-o <binary output file>]
	[--list (list, but do not write binary SRA inventory)]
	[--date.from <YYYY-MM-DD>] (only download SRA records received after this date)
	[--date.to <YYYY-MM-DD>] (only download SRA records received before this date)
	[--strategy <strategy key word>] (only download SRA records that match one of the specified experimental strategies)
		Examples include: RNA-Seq, WGS, AMPLICON, Bisulfite-Seq, ... (case sensitive!)
	[--source <source key word>] (only download SRA records that match one of the specified exterimental sources)
		Examples include: TRANSCRIPTOMIC, GENOMIC, METAGENOMIC, METATRANSCRIPTOMIC, ... (case sensitive!)
	[--include <list of SRA run accessions>] (only download SRA records that match one of the specified SRA runs)
```

- The argument to the `-i` flag is the compressed tar file of XML metadata dowloaded from the NCBI FTP site.
- The argument to the `-o` flag is the name of the binary output file (it is required, unless the the `--list` flag is
  also specified).
- The `--list` flag lists the SRA records found in the XML metafile (and suppressed the creation of the binary output
  file).
- The argument to the `--date.from` flag restricts the download to SRA runs that were *received after* this date
- The argument to the `--date.to` flag restricts the download to SRA runs that were *received before* this date
- The argument to the `--strategy` flag restricts the download to SRA runs that have a matching experimental strategy. 
	- This command can be repeated multiple times.
	- The string matching is case sensitive.
	- A *partial* list of experimental strategies include:
		- RNA-Seq
		- WGS
		- AMPLICON
		- Bisulfite-Seq
- The argument to the `--source` flag restricts the download to SRA runs that have a matching experimental source. 
	- This command can be repeated multiple times.
	- The string matching is case sensitive.
	- A *partial* list of experimental sources include:
		- TRANSCRIPTOMIC
		- GENOMIC
		- METAGENOMIC
		- METATRANSCRIPTOMIC
- The argument to the `--include` flag restricts the download to SRA accession contained within the specified file. The SRA accessions must be listed in an ASCII text file, with one accession per line.

### SRA annotation information

Currently, SRA annotation information is extracted from the XML metadata dowloaded from the NCBI and included in the Bloom filter files and the final CALDERA database. Unless noted, information is stored as a string. Since not all fields are provided for every SRA record, the string "NA" indicates a missing field.
- Run
	- Accession
- Experiment
	- Accession
- Experiment
	- Title
	- Design description
	- Library name
	- Library strategy
	- Library source
	- Library selection
	- Instrument model
- Sample
	- Accession
	- Taxa
	- Attributes (as a collection of key/value pairs)
- Study
	- Accession
	- Title
	- Abstract

## SRA download and database construction

The `maestro` program is responsible for streaming download of SRA records (i.e. downloading SRA files &rightarrow; 
creating Bloom filters &rightarrow; removing SRA files), packing database files with multiple Bloom filters and copying the database files to long-term (e.g. AWS S3) storage. The command line arguments for the `maestro` program are:
```
Usage for Maestro:
	--scratch <parent directory for staging Bloom filter and database files>
	[--scratch.bloom <scratch directory for staging Bloom filter>]
	[--scratch.database <scratch directory for staging database files>]
	--meta <binary metadata input file>
	--s3 <S3 path for storing database files>
	[--s3.no-write (do *not* write database files to s3)]
	[--stream (stream SRA data -- do not use prefetch to download!)]
	[--max-sra-download <max allowed SRA file size in GB>] (default is 30 GB)
	[--status <binary SRA status file for restart>] (default is ./__sra_db_status.bin)
	[--retry <number of download attempts>] (default is 3)
	[--retry.bloom (retry all failed Bloom filters)]
	[--delay <minimum number of seconds between download/streaming requests>]
	[--halt-after <halt after this many SRA downloads> (default is not to stop)]
	[-k <kmer length>] (default is 31)
	[-p <false positive probability (per k-mer, per-filter)>] (default is 0.25)
	[--min-kmer-count <minimum allowed k-mer count>] (default is 5)
	[--hash <hash function name>] (default is murmur32)
		Allowed hash functions: murmur32
	[--len.min <log2 Bloom filter len>] (default is 18)
	[--len.max <log2 Bloom filter len>] (default is 32)
	[-v (turn on verbose output)]
	[--save.bloom (don't remove Bloom filters after database construction)]
	[--save.db (don't remove database file after S3 upload)]
	[--save.sra (don't remove SRA files after Bloom filter construction)]
	[--skip <SRA run accession> (skip over the specified accession; may be repeated)]
```

- The argument to the `--scratch` flag specifies the parent directory for staging intermediate files (e.g. SRA files, Bloom filter files and database files). How much scratch space is needs? Good question! The amount of space is difficult to estimate, as it depends on the number of MPI ranks provided to `maestro`. As a rule of thumb, we recommend providing approximately 2-3 TB of scratch space when running with 150 - 200 MPI ranks. The scratch space must be available to *all* MPI rank and, when running on AWS, the Lustre FSX filesystem appears to work well.
- The argument to the optional `--scratch.bloom` flag specifies the directory for staging intermediate Bloom filter files.
- The argument to the optional `--scratch.database` flag specifies the directory for staging intermediate database files.
- The argument to the `--meta` flag is the binary file of metadata created by the `sra_inventory` program described above. This file defines the SRA accessions that will be downloaded, converted into Bloom filters and stored in database files.
- The argument to the `--s3` flag is the S3 path where the final database files will be stored. Database files are copied to S3 to save expensive local storage and to facilitate the moving of database files to other systems.
- The the optional `--s3.no-write` flag supresses the writing of database files to s3. This avoids the dependancy on AWS command line tools. However, when indexing a large number of SRA records, the space required for storing the database files can be significant, so be prepared!
- The optional `--stream` flag using the NCBI SRA API to stream SRA records (instead of using the NCBI `prefetch` command). While these two approaches are quite similar (both require downloading the "helper" files that are associated with some of the individual SRA records), streamming appears to be slightly more reliable when running on AWS and writing to the FSX/Lustre filesystem. We have observed cases where the `prefetch` tool appears to stall when writing an `.sra` file to FSX/Lustre.
- The argument to the optional `--max-sra-download` flag specifies the maximum size of downloaded SRA records in **GB** (by passing this argument to the `prefetch` command). Limiting the size of SRA downloads is not currently implemented when streaming SRA records using the `--stream` flag.
- The argument to the optional `--status` flag specifies the binary file (by default `./__sra_db_status.bin`) that is used to store the current state of the `maestro` program. This file is **essential** for restarting the `maestro` program after a crash or system failure.
- The argument to the optional `--retry` flag specifies the number of SRA download attempts before giving up on a specific SRA record. The default number of atempts is 3.
- The optional `--retry.bloom` flag forces the program to retry all of the SRA records that previously failed to generate a valid Bloom filter.
- The argument to the optional `--delay` flag specifies the number of seconds to wait between SRA download or streaming requests. By default, there is no delay.
- The argument to the optional `--halt-after` flags specifies the number of SRA records to index (download and convert to a Bloom filter) before halting. This is useful for debugging. By default, the program halts after **all** requested SRA records are indexed.
- The argument to the optional `-k` flag is the integer k-mer length. The default k-mer length is 31 and the maximum allowed k-mer length is 32.
- The argument to the optional `-p` flag is the desired Bloom filter false positive rate (per k-mer and per-filter). The default value is 0.25. This rate must be between 0 and 1, and it used to determine the Bloom filter parameters (i.e. number of bits and number of hash functions used per-k-mer). **SRA files that have more k-mers than can be satisfied by the desired false positive rate threshold (for the specified limits of Bloom filter length and number of hash functions) are skipped**.
- The argument to the optional `--min-kmer-count` flag is the minimum number of times a given k-mer must be found in an SRA file. The default value is 5. The purpose of this threshold is to "de-noise" (i.e. remove sequencing errors from) SRA files. Sequencing errors can generate large numbers of unique, but not biologically relevant k-mers. Larger minimum k-mer counts significantly reduce the size of the Bloom filter needed to represent an SRA file, but also introduce the possibility of false negative search results.
- The argument to the optional `--hash` flag is the name of the hash function to use when constructing the Bloom filter. Currently, only the 32-bit [murmur hash](https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp) function is supported. In the future, additional hash functions may be added (e.g. 64-bit murmur hash).
- The argument to the optional `len.min` flag is the log<sub>2</sub>(minimum Bloom filter length in bits) to test when computing the optimal Bloom filter length. The default value is 18.
- The argument to the optional `len.max` flag is the log<sub>2</sub>(maximum Bloom filter length in bits) to test when computing the optimal Bloom filter length. The defaut value is 32.
- The optional flag `-v` turns on verbose output.
- The optional flag `--save.bloom` prevents the deletion of the intermediate Bloom filter files after database construction. Useful for debugging.
- The optional flag `--save.db` prevents the deletion of the intermediate database files after uploading database files to AWS S3 storage. Useful for debugging.
- The optional flag `--save.sra` prevents the deletion of the intermediate `.sra` files (and associated helper files) after Bloom filter construction. Useful for debugging and only applies when the `--stream` option is **not** specified.
- The argument to the optional `--skip` flag is an ASCII text file of SRA accessions to **skip**. The file must contain one SRA accession per line.

## Searching the database

The `caldera` program enables searching the database files with one or more user-supplied DNA sequences. Currently, all of the database files must be stored on a POSIX filesystem (i.e. **not** on AWS S3).

The command line arguments are:
```
Usage for caldera:
	[-o <output file>] (default is stdout)
	[--o.csv (output CSV) | --o.json (output JSON)]
	[-t <search threshold>] (default is 1)
	-d <database search path> (can be repeated)
	[-i <input sequence file>] (can be repeated)
	[<DNA sequence>] (can be repeated)
```

- The argument to the optional `-o` flag specifies the name of the output file of search results. By default, search results are written to stdout to enable piping to other programs.
	- The output format of the search results can be either comma separated values (`--o.csv`) or JSON format (`--o.json). The default output format is JSON.
- The argument to the optional `-t` flag species the search threshold in terms of the minimum fraction of query sequence kmers that must match a individual database records. Valid threshold values must be greater than 0 and less than, or equal to, 1. The default search threshold is 1 (i.e. **all** query sequence kmers must be found in the same database record to report a match).
- The argument to the `-d` flag specifies the parent directory for all database files. The parent directory and all subdirectories are recursively searched for database files (i.e. `.db` or `.dbz`).
- The argument to the optional `-i` flag is the name of a fasta or fastq file of query sequences that will be searched against all of the specified database files.
- Multiple DNA sequences (delimited by white space) can be directly specified on the command line.
