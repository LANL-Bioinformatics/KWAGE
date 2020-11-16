# Maestro changes
	- Version 0.9 October 19, 2020
		- Improved the buffering of the transposed bitmatrix in build_db() to significantly reduce the number
		  write operations (and hopefully fix the issues with Lustre/FSX stalling/hanging ...).

	- Version 0.8 September 13, 2020
		- Provide for separate scratch paths: `--scratch.bloom` and `--scratch.database`
		- Fixed template specialization-related compiler errors on g++
		- Added the option to `--skip` user-specified SRA run accessions. This is useful when trying to
		  deal with accessions that cannot be downloaded and/or parsed by the SRA toolkit (i.e. ERR1197571).
		  	- At some point, we will need to add an "include" option to be able process previously skipped
			  accessions (should they get fixed and become possible to include again).
		- Fixed bug in final database creation step (when all SRA records have been processed)

	- Version 0.7 August 29, 2020
		- Track worker memory usage to check for memory leaks
		- Even when streaming SRA records, remove *.cache and other SRA record-specific files.

	- Version 0.6 August 15, 2020
		- Handle aligned colorspace SRA records as a special case. These records cannot be read using the
		strategy of loading primary alignments followed by unaligned reads.
		- Use the SRA metadata to extract the number of bases and scale the size of the counting Bloom
		filters (to obtain the desired false positive rate).

	- Version 0.5 August 5, 2020
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

	- Version 0.4 July 30, 2020
		- Modified the restore_bloom() in maestro_main.cpp to restore SRA accessions that are labeled
		  as STATUS_DATABASE_FAIL (in addition to STATUS_BLOOM_SUCCESS). This will recover from previous database creation failures.
		- Moved some file I/O code to a new source file (file_io.cpp) from the maestro_main.cpp file.
		- Created a new help program, "manual_db", to update the accession in a status file with the accessions
		  in a database file that is being manually copied to S3. This is needed when we get an upload failure
		  (due to "aws s3 mv/cp" failing to upload a database file -- most likely due to Lustre FSX being overtaxed
		  during the creation of a large databaes file.)

# SRA inventory changes
	- Version 0.7 September 24, 2020
		- Added the option to specify an optional list of SRA run accessions to include.

# BIGSI++ changes
	- Version 0.4b September 25, 2020
		- Fixed supurious output when the input query is too small for a single kmer.

# Maestro notes and todo
	- DONE: Obtain more information on failed database creation errors; perhaps a separate tool to perform this task ...
	- TODO: Do we need to record the endian-ness in the binary database files?
	
	- Halted database construction to figure out why prefetch is now waiting approximately 4000 sec before returning error 768
		- This appears to be due to an interaction with the Lustre/FSX filesystem provided by AWS
	- Are we obtaining sufficient size reduction? Have the Bloom filter creation code report size reduction. 
		- Use the reported "deflation" value to estimate the size reduction as a funciton of the allowed false positive rate.
		  If total SRA repository is ~5 PB and we would like a 50 TB final database size, we will need a deflation ration of
		  50 TB/5 PB = 50/(5*1024) ~= 0.01
			- The filter length-weighted average deflation is computed using the awk one-liner:
			grep "created Bloom filter" slurm-*.out | awk '{len=2^$13;ave += $22*len;norm += len;++count}END{print "average deflation = " ave/norm "; count = " count}' 
			
			[False positive rate] 	[average deflation]		[sample size]
			0.25					0.0431228				21425
			0.30					0.0613762				64268
			0.35					
			- May need to increase the minimum number of kmers

	- It is difficult to get more than 100 of the m5n.xlarge instance.

	- There are some SRA records that cause the NGS C++ API to leak memory. These may be ABI Solid runs and runs with errors?
		Do we need to try and fix, or are these records rare enough? The slurm-100.log file has %memory consumption information
		that may be useful.
		- ERR2972912

	- It appears that retrying a failed, streaming Bloom filter attempt does not work. A random sampling of the maestro log file
	  shows that almost all (but not every) of the successfull streaming Bloom filter creation steps worked on the first attemp. Why?
	- Is the prefetch timeout issue we're having related to FSX-Lustre? Try pointing the "location of user-respository" to a   *local* directory /scratch/DB --> /tmp
		- Yes -- prefetch seems to have issues with FSX-Lustre. Prefetching files to instance-local storage appears to be faster.
	- Accelerating the Bloom filter calculation:
		- According to the gprof tool, 65% of the calculation time is spent in count_words(), while 10% of the time is spent 
		in murmur_hash32()
	- Try enabling the pagemap thread in SRA/ncbi-vdb/libs/vdb/cursor-table.c
		- ncbi-vdb/libs/vdb/cursor-table.c: The pagemap thread is commented out with an #if 0 -- why?
		- What does ./ncbi-vdb/libs/vdb/cursor-cmn.c: VCursorLaunchPagemapThread() actually do? --> calls launchPagemapThread
		- ncbi-vdb/libs/vdb/cursor-table.c: capacity variable is zero in VTableCreateCachedCursorReadImpl() -- why? (not due to #ifdef)
		- VTableCreateCachedCursorReadImpl() is called via VTableCreateCachedCursorRead() -- however, self->cache_tbl is NULL,
		  so no thread is created.
		  	- VTable *self is the first argument
		- VTableCreateCachedCursorRead() is called by VTableCreateCursorRead() which is called by NGS_CursorMake()
		- ./ncbi-vdb/libs/ngs/NGS_Cursor.c: NGS_CursorMake() -- The table is the *second* argument to this function
			- Called by NGS_CursorMakeDb
				- ./ncbi-vdb/libs/vdb/table-cmn.c : VDatabaseOpenTableRead allocates the table
				- The db variable in NGS_CursorMakeDb(), via self->cache_db, contains the cache_db
		- In ./ncbi-vdb/libs/vdb/database-cmn.c, VDBManagerOpenDBReadVPath() does not set the cache_db variable -- why??
		- KThreadMake ...
	- It appears that streaming API still downloads files. What happens if we use vdb-config to change the filesystem path to 
	use a local directory (instead of Lustre)?

	- Memory leak when streaming SRA data (but not when loading via prefetched-file)
		==346== LEAK SUMMARY:
		==346==    definitely lost: 53,680 bytes in 15 blocks
		==346==    indirectly lost: 8,224 bytes in 2 blocks
		==346==      possibly lost: 0 bytes in 0 blocks
		==346==    still reachable: 1,428,891 bytes in 5,630 blocks
		==346==         suppressed: 0 bytes in 0 blocks
		- In ./libs/kns/http-request.c:
			==346== 4,112 bytes in 1 blocks are definitely lost in loss record 86 of 110
			==346==    at 0x4C2FB0F: malloc (in /usr/lib/valgrind/vgpreload_memcheck-amd64-linux.so)
			==346==    by 0x253B91: KDataBufferResize (in /home/ubuntu/bigsi++3/thread_test)
			==346==    by 0x25ADFD: KDataBufferVPrintf (in /home/ubuntu/bigsi++3/thread_test)
			==346==    by 0x25AF56: KDataBufferPrintf (in /home/ubuntu/bigsi++3/thread_test)
			==346==    by 0x1F0895: KClientHttpRequestFormatMsgInt (in /home/ubuntu/bigsi++3/thread_test)
			==346==    by 0x1F0DB9: KClientHttpRequestSendReceiveNoBody (in /home/ubuntu/bigsi++3/thread_test)
			==346==    by 0x1CEF95: KNSManager_Read (in /home/ubuntu/bigsi++3/thread_test)
			==346==    by 0x1CCD6F: CloudMgrWithinAWS (in /home/ubuntu/bigsi++3/thread_test)
			==346==    by 0x1CBBD0: CloudMgrMake (in /home/ubuntu/bigsi++3/thread_test)
			==346==    by 0x1B8A84: SRequestInitNamesSCgiRequest (in /home/ubuntu/bigsi++3/thread_test)
			==346==    by 0x1B8F75: KService1NameWithVersionAndType (in /home/ubuntu/bigsi++3/thread_test)
			==346==    by 0x1BA325: KService1NameWithVersion (in /home/ubuntu/bigsi++3/thread_test)
		- In KClientHttpRequestSendReceiveNoBodyInt(): 
			- Try adding "KDataBufferWhack( & buffer );" in the max_redirect loop when breaking out of loop and
		  	  when finished with the buffer since it appears that the KDataBuffer is not released in this loop
			- This appears to work!
				==28985== LEAK SUMMARY:
				==28985==    definitely lost: 224 bytes in 2 blocks
				==28985==    indirectly lost: 8,224 bytes in 2 blocks
				==28985==      possibly lost: 0 bytes in 0 blocks
				==28985==    still reachable: 1,428,891 bytes in 5,630 blocks
		- After adding the (apparently) missing KDataBufferWhack() calls to KClientHttpRequestSendReceiveNoBodyInt(), we
		are left with:
			==28985== 4,224 (112 direct, 4,112 indirect) bytes in 1 blocks are definitely lost in loss record 86 of 103
			==28985==    at 0x4C31B25: calloc (in /usr/lib/valgrind/vgpreload_memcheck-amd64-linux.so)
			==28985==    by 0x1EC6FB: KClientHttpAddHeaderString.part.1 (in /home/ubuntu/bigsi++3/thread_test)
			==28985==    by 0x1ED3E9: KClientHttpVAddHeader (in /home/ubuntu/bigsi++3/thread_test)
			==28985==    by 0x1ED47C: KClientHttpAddHeader (in /home/ubuntu/bigsi++3/thread_test)
			==28985==    by 0x1F30A1: KClientHttpRequestPOST_Int (in /home/ubuntu/bigsi++3/thread_test)
			==28985==    by 0x1F3374: KClientHttpRequestPOST (in /home/ubuntu/bigsi++3/thread_test)
			==28985==    by 0x1F3574: KClientHttpRequestHEAD (in /home/ubuntu/bigsi++3/thread_test)
			==28985==    by 0x1EB811: KNSManagerVMakeHttpFileIntUnstableImpl (in /home/ubuntu/bigsi++3/thread_test)
			==28985==    by 0x1EC560: KNSManagerVMakeHttpFileIntUnstableFromBuffer (in /home/ubuntu/bigsi++3/thread_test)
			==28985==    by 0x1E5B06: KNSManagerMakeReliableHttpFile (in /home/ubuntu/bigsi++3/thread_test)
			==28985==    by 0x1AC46F: VFSManagerMakeHTTPFile (in /home/ubuntu/bigsi++3/thread_test)
			==28985==    by 0x1AECF9: VFSManagerOpenDirectoryReadDecryptRemote (in /home/ubuntu/bigsi++3/thread_test)
	
	-Question for SRA -- can we know in advance how big the SRA cache files will be? There can be large
		cache files, for example:
		-rw-rw-r-- 1 ubuntu ubuntu  34G Sep 11 14:14 SRR628885.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  36G Sep 11 14:14 SRR628887.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  34G Sep 11 14:14 SRR628888.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  61G Sep 12 07:30 SRR650708.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  44G Sep 12 07:30 SRR651964.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  41G Sep 12 07:30 SRR652436.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  31G Sep 12 07:30 SRR654044.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  30G Sep 12 07:30 SRR654045.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  31G Sep 12 07:30 SRR654046.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  32G Sep 12 07:30 SRR654047.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  33G Sep 12 07:30 SRR663248.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  31G Sep 12 07:30 SRR663452.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  22G Sep 12 07:30 SRR671683.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  25G Sep 12 07:30 SRR674285.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  26G Sep 12 07:30 SRR674482.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  24G Sep 12 07:30 SRR675160.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  13G Sep 12 07:30 SRR675233.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  24G Sep 12 07:30 SRR680161.sra.cache
		-rw-rw-r-- 1 ubuntu ubuntu  22G Sep 12 07:30 SRR680163.sra.cache
	- How to change the location of the SRA cache files independent of the RefSeq directory?
	- Amazon appears to charge $0.1 per GB per month - how much of a benefit do we get from storing SRA cache files locally?

		Costs for first 12 days in September:
			Amazon Elastic Compute Cloud NatGateway	 = $61.90
				$0.045 per GB Data Processed by NAT Gateways
				815.445 GB
				$36.70
				$0.045 per NAT Gateway Hour
				560.000 Hrs
				$25.20
			Amazon Elastic Compute Cloud running Linux/UNIX = $4,493.47
				$0.085 per On Demand Linux c5.large Instance Hour
				280.000 Hrs
				$23.80
				$0.126 per On Demand Linux r5.large Instance Hour
				34,513.580 Hrs
				$4,348.71
				$0.432 per On Demand Linux c5n.2xlarge Instance Hour
				280.000 Hrs
				$120.96
			EBS = $170.68
				$0.10 per GB-month of General Purpose SSD (gp2) provisioned storage - US East (Northern Virginia)
				1,706.844 GB-Mo
				$170.68
	- What if we allocate 500 GB of EBS storage per node and write all SRA data (including refseq) and database files (prior to
	  sending to an S3 bucket) locally?
		- Can we stream the writing of database files to S3? We could then get away with approx 200 GB of local storage.
		- Limit the number of Bloom filters per database file for larger Bloom filters. We had used 2048 filters/per file
		  to enable efficient bitslice compression
		  	- Since a 2048 bit slice requires 256 uncompressed bytes to store, the number of compressed bytes ( <= 256) can be
			  stored in a single byte.
		- log L = 32 -> 2^32 filter bits
		- There are 2048 filters/database and each bitslice requires 256 bytes (i.e. 2048/8)
		- Total number of bytes to store a datbase file = 2^32 * 256 = 1024 GB = 1 TB <-- Yiikes! Too big
			- according to S3, there are some 256 GB files: 2020-09-06 23:24:01  256.0 GiB sra.431.db
		- In addition to specifying the number of filters per file *also* specify the maximum database file size
	
	- The following records appears to be "pathological" in that they crash nodes that try to read them and 
	  compute Bloom filters (consuming all the RAM?). Are these all aligned colorspace reads (which the SRA-team have
	  already said are "broken" an unlikely to be fixed)?
	  	SRR579581
	  	ERR661170
	  	SRR643883
		SRR819819
		SRR891318
		ERR1197571 <- Causes rapid crash

			prefetch ERR1197571                                                                  

			2020-09-26T17:08:38 prefetch.2.10.8: 1) Downloading 'ERR1197571'...
			2020-09-26T17:08:38 prefetch.2.10.8:  Downloading via HTTPS...
			2020-09-26T17:08:45 prefetch.2.10.8:  HTTPS download succeed
			2020-09-26T17:08:46 prefetch.2.10.8 int: no error - failed to verify
			2020-09-26T17:08:46 prefetch.2.10.8: 1) failed to download ERR1197571

		SRR8388771 <- Causes rapid crash

			prefetch SRR8388771

			2020-10-02T20:42:24 prefetch.2.10.8: 1) Downloading 'SRR8388771'...
			2020-10-02T20:42:24 prefetch.2.10.8:  Downloading via HTTPS...
			2020-10-02T20:42:26 prefetch.2.10.8:  HTTPS download succeed
			2020-10-02T20:42:27 prefetch.2.10.8 int: no error - failed to verify
			2020-10-02T20:42:27 prefetch.2.10.8: 1) failed to download SRR8388771
	
	- Need to add the ability to perfom global compression of database files.
		- Rather than compress individual bitslices (as would be needed for opertational compression), we can compress
		  entire database files to reduce the cost of transmission and long term storage. For example:
		  Before compression: -rw-rw-r-- 1 ubuntu ubuntu 4.1G Sep 25 00:34 sra.1.db
		  After comprresion : -rw-rw-r-- 1 ubuntu ubuntu 2.4G Sep 25 00:34 sra.1.db.gz
		
		- This is likely only to help with the prioritized SRA data (which has many SRA files for a small number
		  of taxonomic groups).
		- For now, write a script to copy S3 objects to local disk -> compress -> copy compressed file back to S3 
		  (as a separate object to minimize risk of losing any database files).

	- Priority list of SRA database files is 1.8 TB (uncompressed) in 744 files
		- After merging files: 1.8 TB 591 files
		- After compressing: 1.1 TB 591 files
# Name change from BIGSI++ to
	- Sriracha
	- Falchion (lots of github projects)
	- Curtana (lot of github projects)
	- Poignard / Poniard
	- ecilstib ("bitslice" backwards)
	- Sawfish (a window manager ...)
	- Xiphias (genus of sword fish, some small github projects)
	- Elasmobranch (class containing sharks)
		- Elasmosearch (taken!)
	- Seqflower
	- shovelnose (the genus that contains the Sawfish)
	- Bloom Based Bit Sliced Sequence Search Algorithm
	- zBLASE: compressed "complressed BLoom based Algorithm for Sequence Explorations"
	- PERFUME: ???
	- indexium (German company)
	- key words:
		- Bloom
		- Binary, Bit
		- Sequence
		- Search
		- Filter
		- Parallel
		- Slice, Cut, Chop
	- Parts of a flower:
		- pistil
		- stamen
		- stigmata
		- anther
		- filament
	- BIGSI == "big sigh"
		- RELAXED
		- zBLISS <- Already used by a bioinformatics tool ...
			- compressed BLoom fIlter Sequence Search
		- 


	
