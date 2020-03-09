# BIGSI++
A C++ reimplementation of the ultra-fast, [Bloom filter](https://en.wikipedia.org/wiki/Bloom_filter)-based sequence search originally developed by the [Iqbal group](https://www.nature.com/articles/s41587-018-0010-1). This project is similar in spirit to the Iqbal group's BIGSI follow-on project: [COBS: a Compact Bit-Sliced Signature Index](https://arxiv.org/abs/1905.09624). While the goals are very similar, the BIGSI++ project is an independent implementation of the original BIGSI algorithm with the addition of the following features:
1. Parallel, C++ implementation (using both OpenMP and MPI)
2. Adaptive Bloom filter size (similar to COBS)
3. Compressed Bloom filters
4. Bloom filter construction directly from [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) (SRA) files

The BIGSI++ pipeline has four primary components:
1. Streaming downloads of SRA files
2. Streaming construction of Bloom filters from the downloaded SRA files
3. Database construction from the Bloom filters
4. Searching the resulting database with a nucleic acid query

## Building and installing the code
The BIGSI++ code is written in C++ and requires the MPI (message passing interface; I recommend [OpenMPI](https://www.open-mpi.org/)). A compiler that supports OpenMP is also suggested, but not required.

Currently, downloading the SRA files is orchestrated by the `sra_download` program, with help from the `prefetch` tool (that is provided in the SRA toolkit) and a cluster scheduling system (so that Bloom filters can be created from SRA files while new files are being downloaded in parallel). The code currently assumes a [Torque/PBS](https://www.mcs.anl.gov/research/projects/openpbs/) compatible scheduler, but this is only for the current proof-of-principle implementation. As discussed below, downloading the entire SRA to a private network (i.e. non-[AWS](https://aws.amazon.com/), non-[Google](https://cloud.google.com/)) is **not** going to work -- downloading will take too long for the several petabytes of SRA data.

A cloud-compatible version of `sra_download` (coming soon!) will replace the Torque/PBS dependency with an AWS dependency to spawn additional AWS instances than compute Bloom filters from the downloaded files.

### Download the SRA toolkit and associated APIs
In addition, the SRA toolkit and associated APIs are required if you would like to be able to download and/or construct Bloom filters from SRA files. In particular, you will need to download and build the following three GitHub projects from NCBI:
1. [ncbi-vdb](https://github.com/ncbi/ncbi-vdb)
2. [ngs](https://github.com/ncbi/ngs)
3. [sra-tools](https://github.com/ncbi/sra-tools)

Please note that these NCBI libraries are only needed if you will need to either download SRA files and/or create Bloom filters from SRA files. If you only need to either build a database from existing Bloom filters or search an existing database, then these NCBI libraries are *not* needed.

### Edit the BIGSI++ Makefile
After the three separate NCBI software packages have been downloaded and compiled, you will need to edit the BIGSI++ Makefile to specify the locations of both the include and library directories for the ncbi-vdb and ngs packages.

Edit the `SRA_LIB_PATH` variable to point to the directory on your system that contains the SRA library files.

Edit the `SRA_INCLUDE_PATH` variable to point to the directory on your system that contains the SRA include files.

## SRA download

Streaming download of SRA files (i.e. download SRA file &rightarrow; create Bloom filter &rightarrow; remove SRA file) is performed by the `sra_download` program, whose command line arguments are:
```
Usage for SRA Download:
	-i <XML metadata file from ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata>
	--download <download directory for SRA data>
	--bloom <Bloom filter output directory>
	--log <log filename> (required to track progress for restarts)
	[-k <kmer length>] (default is 31)
	[-p <false positive probability (per k-mer, per-filter)>] (default is 0.25)
	[--min-kmer-count <minimum allowed k-mer count>] (default is 5)
	[--hash <hash function name>] (default is murmur)
		Allowed hash functions: murmur
	[--list (list, but do not download, SRA data)]
	[--sleep <sec> (time to sleep between downloads)]
	[--max-backlog <number of downloads> (Pause SRA downloading when exceeded)] (default is 25; 0 is no limit)
	[-t <number of download threads>] (default is 1)
```

- The argument to the `-i` flag is the compressed tar file of XML metadata dowloaded from the NCBI [FTP](tp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata) site. **Do not uncompress or untar this file!** Since the uncompressed/untarred files are large and cumbersome to manage, the `sra_download` program directly reads the SRA metadata and run accessions by uncompressing on the fly (and parses the tar file format to track the source files).
- The argument to the `--download` flag is the directory that will be used to temporarily stage SRA files during the download process and while the Bloom filters are being created by the `bloomer` tool. Ideally, this directory should be located on a fast storage device and must be visible to all of the backend node/cloud instances that will handle the conversion of SRA files into Bloom filters.
- The argument to the `--bloom` flag is the directory that will be used to store all of the Bloom filters.
- The argument of the `--log` flag is the filename of either a new log file to create (for a completely new run) or an existing log file (for restarting a download attempt or incrementally downloading newly available SRA files). This log file is a human readable text file that is also parsed by the `sra_download` program to identify which SRA files have already been downloaded.
- The argument to the `-k` flag is the integer k-mer length. The maximum allowed k-mer length is 32.
- The argument to the `-p` flag is the desired Bloom filter false positive rate (per k-mer and pe-filter). This rate must be between 0 and 1, and it used to determine the Bloom filter parameters (i.e. number of bits and number of hash functions used per-k-mer). Currently, the maximum allowed Bloom filter length is 2<sup>31</sup> bits (which makes the maximum file size of a single Bloom filter approximately 256 MB). SRA files that have so many k-mers that the desired false positive rate threshold can not be obtained without exceeding this maximum Bloom filter length are skipped (and an empty Bloom filter file is written). In the future, the maximum Bloom filter size will be a user configurable option (and will support a larger number of bits).
- The argument to the `--min-kmer-count` flag is the minimum number of times a given k-mer must be found in an SRA file. The purpose of this threshold is to "de-noise" (i.e. remove sequencing errors from) SRA files. Since sequencing errors can generate larg numbers of unique, but not biologically relevant k-mers, larger minimum k-mer counts significantly reduce the size of the Bloom filter needed to represent an SRA file.
- The argument to the `--hash` flag is the name of the hash function to use when constructing the Bloom filter. Currently, only the [murmur hash](https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp) function is supported. In the future, additional hash functions may be added. I suspect that there is still room to improve hash function performance (i.e. increase speed, reduce collisions) for DNA sequences -- however, hash function construction is very tricky!
- Specifing the `--list` flag will cause `sra_download` to print the SRA run accessions and metadata found in the XML input file and then exit **without** downloading any SRA files.
- The argument to the optional `--sleep` flag is the number of seconds to sleep after downloading an SRA file using the `prefetch` command. This is intended for situations where one might need to reduce the average bandwidth to NCBI SRA database.
- The argument to the `--max-backlog` flag sets the maximum number of SRA files that can be downloaded before being processed with the `bloomer` program to create Bloom filters (and to remove the SRA file). This is intended to prevent the exhaustion of local `--download` directory storage in the event that (a) the job scheduler has failed and/or (b) the number of compute nodes running the `bloomer` program can not keep up with the rate at which SRA files are being downloaded.
- The argument to the `-t` command sets the integer number of threads that will be used to download SRA files in parallel (i.e. spawning independent instances of the `prefetch` program to download separate SRA run accessions). When downloading to a private (i.e. non-cloud) network, using upto 4 threads offers a bandwidth improvement.


## Bloom filter construction

## Database construction

## Searching the database
