import pysam
import sys, getopt
import re
import commands

def getFIlesInuput(argv):
	output={'ChimericFile': '', 'out': ''}
	try:
		opts, args = getopt.getopt(argv, "i:c:o:")
	except getopt.GetoptError:
		sys.exit(2)
	for opt, arg in opts:
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
	file_Chimeric = open(myFiles['ChimericFile'], 'r')
	file_out = open(myFiles['out'], 'w')
	for line in file_Chimeric.readlines():
		colmunsLine = re.split(r'\t', line)
		if colmunsLine[0] != colmunsLine[3] or colmunsLine[2] != colmunsLine [5]:
			file_out.write("\t".join(colmunsLine))
		else :
			StartSegment, StartCIGAR, EndSegment, EndCIGAR = (colmunsLine[10], colmunsLine[11], colmunsLine[12], colmunsLine[13]) if colmunsLine[2] == "+" else (colmunsLine[12], colmunsLine[13], colmunsLine[10], colmunsLine[11])
			shiftAtEnd = findNumberOFShift(EndCIGAR)
			if (int(EndSegment) + shiftAtEnd) < int(StartSegment) :
				file_out.write("\t".join(colmunsLine))

			else :
				continue
	file_Chimeric.close()
	file_out.close()

main()



