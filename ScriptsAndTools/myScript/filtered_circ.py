#!/usr/bin/python2.7

import sys, getopt
import re

def getFIlesInuput(argv):
	output={'splice_fa': '', 'splice_bed': '', 'path_out': ''}
	try:
		opts, args = getopt.getopt(argv, "b:f:p:")
	except getopt.GetoptError:
		print("filter -f <splice.fa>")
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-f" :
			output['splice_fa'] = arg
		if opt == "-b" :
			output['splice_bed'] = arg
		if opt == "-p" :
			output['path_out'] = arg
	return output;

def main():
	myFiles = getFIlesInuput(sys.argv[1:])
	file_fasta =  open(myFiles['splice_fa'], "r")
	file_bed = open(myFiles['splice_bed'], "r")
	filtered_fasta = open(myFiles['path_out']+"_filtered_splice_circ.fa", "w")
	filtered_bed = open(myFiles['path_out']+"_filtered_splice_circ.bed", "w")
	boolean=1
	for line in file_fasta.readlines():
		if re.match(r'^>', line) and "_circ_" in line:
			boolean=1
		elif "_norm_" in line:
			boolean=0
		if boolean == 1:
			filtered_fasta.write(line)	


	file_fasta.close()
	filtered_fasta.close()
	for line in file_bed.readlines():
		if "CIRCULAR" in line:
			filtered_bed.write(line)

	file_bed.close()
	filtered_bed.close()
main()


