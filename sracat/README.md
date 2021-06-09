# sracat (version 0.1, June 9, 2021): a command-line tool for extracting unordered read data from SRA files

`sracat` is small C++ program that uses the [NCBI sra-toolkit](https://github.com/ncbi/sra-tools) C++ API to extract sequence (and optionally quailty scores) from SRA records.
Unlike the `fasterq-dump` program included with the sra-toolkit, sracat does *not* output the reads in their origianlly submitted order, but rather outputs
reads in the order in which they are stored in the SRA file.


## Usage

The `sracat` program has the following command-line options:

- `[-o <output file name>]` (by default, output is written to stdout)
- `[--qual]` (write both sequence and quality scores in **fastq** format. By default, only sequence data is written in **fasta** format)
- `<SRA accession or filename 1> ...` (any number of SRA accessions or filenames may be specified -- all output will be concatinated)

`sracat` will read local SRA files that have a `.sra` file extension.

# Caveats
- For SRA records that contain read pairs aligned to a reference, `sracat` will only output read pairs that are either both aligned to the reference or both unaligned to the reference. Pairs with one aligned and one unaligned read will **not** be output!
- For SRA records that contain read pairs aligned to a reference, `sracat` will only generate fasta output (not fastq). The `--qual` command line option will be ignored for these records.
