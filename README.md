# BIGSI++
A C++ reimplementation of the ultra-fast, Bloom filter-based sequence search originally developed by the [Iqbal group](https://www.nature.com/articles/s41587-018-0010-1). This project is similar in spirit to the Iqbal group's BIGSI follow-on project: [COBS: a Compact Bit-Sliced Signature Index](https://arxiv.org/abs/1905.09624). While the goals are very similar, the BIGSI++ project is an independent implementation of the original BIGSI algorithm with the addition of the following features:
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

### Download the SRA toolkit and associated APIs
In addition, the SRA toolkit and associated APIs are required to download and read SRA files. In particular, you will need to download and build the following three GitHub projects from NCBI:
1. [ncbi-vdb](https://github.com/ncbi/ncbi-vdb)
2. [ngs](https://github.com/ncbi/ngs)
3. [sra-tools](https://github.com/ncbi/sra-tools)

### Edit the BIGSI++ Makefile
After these three separate NCBI software packages have been downloaded an compiled, you will need to edit the BIGSI++ Makefile to specify the locations of both the include and library directories for the ncbi-vdb and ngs packages.

## SRA download

## Bloom filter construction

## Database construction

## Searching the database
