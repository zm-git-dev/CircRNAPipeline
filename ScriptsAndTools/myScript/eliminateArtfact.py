import pysam
import sys, getopt
import re
import commands

def getFIlesInuput(argv):
	output={'mergeFile': '', 'ChimericFile': '', 'out': ''}
	try:
		opts, args = getopt.getopt(argv, "i:c:o:")
	except getopt.GetoptError:
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-i" :
			output['mergeFile'] = arg
		if opt == "-c" :
			output['ChimericFile'] = arg
		if opt == "-o" :
			output['out'] = arg
	return output;





def findNumberOFShift(cigar):

	decomposedCIGAR = list(re.findall('(\d+)([MIDNSHPp=XB])', cigar))
	#print(decomposedCIGAR)
	shiftCount = 0
	indiceCIGAR = "".join([x[1] for x in decomposedCIGAR])
	valueCIGAR = [x[0] for x in decomposedCIGAR]
	# On ne recupere les valeurs du N seulement si on a un autre match apres 
	for m in re.finditer('([Np])', indiceCIGAR) :
		if re.search('M', indiceCIGAR[m.start():]) :
			shiftCount += int(valueCIGAR[m.start()])
	#On recupere toute les valeur de M
	for m in re.finditer('M', indiceCIGAR) :
		shiftCount += int(valueCIGAR[m.start()])
	
	return shiftCount


			
def main():
	myFiles = getFIlesInuput(sys.argv[1:])
	file_Merge = open(myFiles['mergeFile'])
	# file_Chimeric = open(myFiles['ChimericFile'])
	file_out = open(myFiles['out'], 'w')
	for line in file_Merge.readlines():
		colmunsLine = re.split(r'\t', line)
		columnCircRnaF = set(re.split(r',', colmunsLine[9])[:-1])
		columnCircExplo = set(re.split(r',', colmunsLine[11])[:-1])
		#print(colmunsLine)
		reads = list(columnCircRnaF | columnCircExplo)
		if len(reads) == 0 :
			file_out.write("\t".join(colmunsLine))
		else :
			listOfCorrectReads = list()
			for read in reads :
				cmd = "grep '" + read + "' " + myFiles['ChimericFile']
				test = commands.getoutput(cmd)
				cmdResult = list(re.split(r'\t', test))
				StartSegment, StartCIGAR, EndSegment, EndCIGAR = (cmdResult[10], cmdResult[11], cmdResult[12], cmdResult[13]) if cmdResult[2] == "+" else (cmdResult[12], cmdResult[13], cmdResult[10], cmdResult[11])
				# print(cmdResult);
				#print(cmdResult[2])
				#print(StartSegment + "\t" +  StartCIGAR + "\t" + EndSegment + "\t" + EndCIGAR)						
				shiftAtEnd = findNumberOFShift(EndCIGAR)
				if (int(EndSegment) + shiftAtEnd) < int(StartSegment) :
					listOfCorrectReads.append(read)

			if len(listOfCorrectReads) > 0 :			
				newColumnCircRNAF = list()
				newColumnCircExplo = list()
				for correctRead in listOfCorrectReads :
					if correctRead in columnCircRnaF :
						newColumnCircRNAF.append(correctRead)
					if correctRead in columnCircExplo :
						newColumnCircExplo.append(correctRead)

				colmunsLine[9] = ",".join(newColumnCircRNAF)
				colmunsLine[11] = ",".join(newColumnCircExplo)
				colmunsLine[8] = str(len(newColumnCircRNAF))
				colmunsLine[10] = str(len(newColumnCircExplo))
				#print(colmunsLine)
				file_out.write("\t".join(colmunsLine))	

	file_Merge.close()
	file_out.close()

main()



