#!/usr/bin/python2.7

import sys, getopt
import re
import commands

def getFIlesInuput(argv):
	output={'LogOutStar': '', 'STARAlignedBam': '', 'Bowtie2AlignedBam': '', 'anchorsFile':'', 'outputPath':'', 'filterChimeric':''}
	try:
		opts, args = getopt.getopt(argv, "l:s:b:a:o:c:")
	except getopt.GetoptError:
		print("LogOutStar -l <pathToLogFinalStar>\nSTARAlignedBam -s <pathToSTARAlignedBam>\nBowtie2AlignedBam -b <pathToBowtie2AlignedBam>\n-a anchorsFile <pathToAnchorsFile>\n -o outputPath <pathTooutput>\n -c <pathToFilterChimericReads>\n")
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-l" :
			output['LogOutStar'] = arg
		if opt == "-s" :
			output['STARAlignedBam'] = arg
		if opt == "-b" :
			output['Bowtie2AlignedBam'] = arg
		if opt == "-a" :
			output['anchorsFile'] = arg
		if opt == "-o" :
			output['outputPath'] = arg
		if opt == "-c" :
			output['filterChimeric'] = arg
	return output;

def main():
	Infos = getFIlesInuput(sys.argv[1:])
	# Open File :
	fileOutput = open(Infos['outputPath'], "w")

	# Get NbReads at Inputs : 
	cmdNbReads = "grep 'input reads' " + Infos['LogOutStar'] +  " | cut -d\| -f2"
	NbReadsLine = commands.getoutput(cmdNbReads)
	mNbReads = re.search('\s+(\d+)*', NbReadsLine)
	NbReads = mNbReads.group(1)
	print(NbReads)
	# Get NbReads Mapped STAR
	cmdNbReadsMappedSTAR = "samtools view -F4 "+ Infos['STARAlignedBam'] + " | cut -f 1 | sort | uniq | wc -l "
	# print(cmdNbReadsMappedSTAR)
	NbReadsMappedSTAR = commands.getoutput(cmdNbReadsMappedSTAR)
	print(NbReadsMappedSTAR)
	
	# Get NbReads Mapped Bowtie2
	cmdNbReadsMappedBowtie2 = "samtools view -F4 " + Infos['Bowtie2AlignedBam']  + " | cut -f 1 | sort | uniq | wc -l "
	# print(cmdNbReadsMappedBowtie2)
	NbReadsMappedBowtie2 = commands.getoutput(cmdNbReadsMappedBowtie2)
	print(NbReadsMappedBowtie2)

	# Get Number Of Spliced Reads STAR
	cmdSplicedReadsSTAR = "samtools view " + Infos['STARAlignedBam']  + " | awk '($6 ~ /N/)' | wc -l "
	# print(cmdSplicedReadsSTAR)
	NbSplicedReadsStar = commands.getoutput(cmdSplicedReadsSTAR)
	print(NbSplicedReadsStar)

	# Get Number Of Spliced Reads bowtie2
	cmdSplicedReadsBowtie2 = "samtools view " + Infos['Bowtie2AlignedBam']  + " | awk '($6 ~ /N/)' | wc -l "
	# print(cmdSplicedReadsBowtie2)
	NbSplicedReadsBowtie2 = commands.getoutput(cmdSplicedReadsBowtie2)
	print(NbSplicedReadsBowtie2)
	
	# Get Number of chimeric reads
	cmdNbChimericReads = "grep 'Number of chimeric reads' " + Infos['LogOutStar']  + " | cut -d\| -f2"
	# print(cmdNbChimericReads)
	NbChimericReadsLine = commands.getoutput(cmdNbChimericReads)
	mNbChimericReads = re.search('\s+(\d+)*', NbChimericReadsLine)
	NbChimericReads = mNbChimericReads.group(1)
	print(NbChimericReads)
	
	# Get Number of Anchors
	cmdNbAnchors = "zcat " + Infos['anchorsFile'] + " | wc -l"
	# print(cmdNbAnchors)
	NbAnchors = int(commands.getoutput(cmdNbAnchors))
	NbAnchors = NbAnchors/4
	print(NbAnchors)	
	
	# Get Number Of Filtered Chimerics Reads
	cmdFilteredChimericReads = "wc -l " + Infos['filterChimeric'] + " | awk '{print $1}' "
	NbFilteredChimeric = commands.getoutput(cmdFilteredChimericReads)
	print(NbFilteredChimeric)

	# Get Number Of Unique Chemeric junctions
	cmdUniqueChimericJunctions = "awk '{print $1$2$3$4$5$6}' " + Infos['filterChimeric'] +  " | sort | uniq | wc -l "
	NbUniqueJunctions = commands.getoutput(cmdUniqueChimericJunctions)
	print(NbUniqueJunctions)
	
	# Print file
	fileOutput.write("NbReadsTotal\tNbReadsMappedBySTAR\tNbReadsMappedByBowtie2\tNbSplicedReadsSTAR\tNBSplicedReadsBowtie2\tNbChimericReads\tNbAnchors\tNbFilteredChimericReads\tNbUniqueChimericJunctions\n")
	fileOutput.write(NbReads + "\t" + NbReadsMappedSTAR + "\t" + NbReadsMappedBowtie2 + "\t" + NbSplicedReadsStar + "\t" + NbSplicedReadsBowtie2 + "\t" + NbChimericReads + "\t" + str(NbAnchors) + "\t" + NbFilteredChimeric + "\t" + NbUniqueJunctions + "\n")
	fileOutput.close()

main()
	
