#!/usr/bin/env python3
# Linux: /home/jgans/anaconda3/bin/python3
# OS X: /anaconda3/bin/python3
#
# Version 0.1 (May 12, 2020):
#	- Initial implementation

import sys
import os
import re
import time
import subprocess

class Record:
	def __init__(self, accession = "", month = 0, library_strategy = "", library_source = ""):
                
		self.accession = accession
		self.month = month
		self.library_strategy = library_strategy
		self.library_source = library_source
	
	def __eq__(self, rhs):

		return self.accession == rhs.accession

def main():

	output_dir = "output"
	search_subject = "NC_045512.fna"

	# Slurm parameters
	num_nodes = 2
	task_per_node = 8
	cpus_per_task = 1

	# The number of jobs to allow in the Slurm queue (both working and running)
	queue_size = 64

	if ( len(sys.argv) != 2 ):

		sys.stderr.write( "Usage: {} <csv inventory file>\n".format(sys.argv[0]) )
		sys.exit(0)

	# Parse the input CSV file
	fin = open(sys.argv[1])

	if not fin:

		sys.stderr.write( "Unable to open CSV inventory file\n".format(sys.argv[1]) )
		sys.exit(0)

	sra = list()

	for line in fin:

		data = line.split(',')

		# No header in the CSV, but we expect the following columns:
		#	accession,date,strategy,source,bases
		if len(data) != 5:

			sys.stderr.write( "Did not find the expected number of columns: {}\n".format(line) )
			sys.exit(0)

		accession = data[0]
		date = data[1]

		m = re.search('\d+-(\d+)-\d+', date)

		month = int( m.group(1) )

		if (month < 1) or (month > 12):
			
			sys.stderr.write("Month is out of bounds\n")
			sys.exit(0)

		strategy = data[2].upper().replace(' ', '_')
		source = data[3].upper().replace(' ', '_')

		if strategy == '':
			strategy = 'UNKNOWN'
		
		if source == '':
			source = 'UNKNOWN'

		sra.append( Record(accession, month, strategy, source) )
	
	fin.close()

	num_accession = len(sra)

	print("Found a total of", num_accession, "records")

	# Count the number of different strategies and sources
	strategy = dict()
	source = dict()

	working = list()
	complete = list()

	for r in sra:

		if(r.library_strategy not in strategy):
			strategy[r.library_strategy] = 0
		
		if(r.library_source not in source):
			source[r.library_source] = 0

		strategy[r.library_strategy] = strategy[r.library_strategy] + 1
		source[r.library_source] = source[r.library_source] + 1

	print("Library strategies:")
	print_sorted_dict(strategy)

	print("Library sources:")
	print_sorted_dict(source)
	
	# Submit each record to the slurm scheduler
	# Write the output files to a directory tree
	# organized by:
	#	output_dir/
	#		SOURCE/
	#			STRATEGY/
	#				MONTH/
	#					accession1
	#					accession2
	#					...

	if not os.path.isdir(output_dir):
		os.mkdir(output_dir)

	num_skipped = 0

	while len(sra) > 0:

		r = sra.pop()

		# Make sure we have a valid path
		path = output_dir + "/" + r.library_source

		if not os.path.isdir(path):
			os.mkdir(path)

		path = path + "/" + r.library_strategy

		if not os.path.isdir(path):
			os.mkdir(path)
		
		path = path + "/" + str(r.month)

		if not os.path.isdir(path):
			os.mkdir(path)

		# Is there an accession file?
		filename = path + "/" + r.accession + ".txt"

		if search_complete(filename):

			num_skipped = num_skipped + 1
			continue

		# DEBUG
		print("Submitting", r.accession)

		script = "#!/bin/bash\n" \
			"#SBATCH --nodes={}\n" \
			"#SBATCH --ntasks-per-node={}\n" \
			"#SBATCH --cpus-per-task={}\n" \
			"#SBATCH --output=/dev/null\n" \
			"\n" \
			"cd /home/ubuntu/jgans/sriracha\n" \
			"\n" \
			"mpirun ./sriracha -i {} \\\n" \
			"		-o {} \\\n" \
			"		--retry 3 \\\n" \
			"		-t 0.5 -n 25 {}\n".format(num_nodes, task_per_node, cpus_per_task, search_subject, filename, r.accession)

		sbatch = subprocess.Popen('/Users/jgans/src/BIGSI/stream/sbatch', stdin=subprocess.PIPE)

		sbatch.stdin.write( script.encode() )
		sbatch.communicate()

		if sbatch.returncode != 0:

			# When the filesystem runs out of space, sbatch cannot write the script to disk and fails.
			# While there may be some other non-fatal failure mode, I haven't found it (so just give up)
			sys.stderr.write( "sbatch encountered an error submitting {}\n".format(r.accession) )
			sys.exit(0)

		working.append(r)

		# Debug
		#print("Pushed", r.accession, "to working")

		if len(working) < queue_size:
			continue

		# Wait until at least one of the working jobs is done
		sleep_count = 0

		while True:

			reaper = list()

			for r in working:

				filename = output_dir + "/" + r.library_source + "/" + r.library_strategy + "/" + str(r.month) + "/" + r.accession + ".txt"

				if not os.path.isfile(filename):
					continue
			
				if search_complete(filename):
					reaper.append(r)
			
			for r in reaper:

				working.remove(r)
				complete.append(r)

				if sleep_count > 0:

					print()

					# Increment the sleep count to 1 to stop printing blank lines
					# for every reaped genome
					sleep_count = 1

				print(r.accession, " done; ", (100.0*len(complete))/(num_accession - num_skipped),"% complete", sep='')

			if len(reaper) == 0:

				if sleep_count == 0:
					print("Sleeping",end='', flush=True)
				else:
					print(".", end='', flush=True)

				time.sleep(10)

				sleep_count = sleep_count + 1
			else:
				break

		#sys.stderr.write( "DEBUG\n" )
		#sys.exit(0)

		# Schedule the search by running the sbatch command and
		# piping the script to stdin

def search_complete(m_filename):

	if not os.path.isfile(m_filename):
		return False

	fin = open(m_filename)

	if not fin:
		sys.stderr.write( "search_complete: Unable to open file {}\n".format(m_filename) )
		sys.exit(0)

	for line in fin:

		if '//' in line:
			return True

	fin.close()

	return False

def print_sorted_dict(m_dict):

	_tmp = list()

	for k,v in m_dict.items():
		_tmp.append( (v, k) )
	
	_tmp.sort()
	_tmp.reverse();

	for x in _tmp:
		print('\t', x[0], '\t', x[1], sep='')

main()