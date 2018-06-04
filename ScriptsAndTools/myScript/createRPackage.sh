#!/bin/bash


# g : genome : path to the genome file
# r : respository : path to the final repository where the package will be created

# How to USE : BSgenome.Hg38.UCSC.1

while getopts ":g:r:s:o:n:" option; do
	case "${option}" in 
	g)
		genome=${OPTARG}
		;;
	r)
		outputRepository=${OPTARG}
		;;
	o)
		organism=${OPTARG}
		;;
	s)
		source=${OPTARG}
		;;
	n)
		namesSeq=${OPTARG}
		;;
	esac
done

GenomeFilename=$(basename $genome)
GenomeName="${GenomeFilename%.*}"

echo $GenomeName


mkdir -p $outputRepository
awk -v outputR=$outputRepository '{if ($1 ~ /^>/){print $1; split($1, seqname, ">"); print seqname[2]; printf $0"\n" >> outputR "/" seqname[2] ".fa"} else{printf $0"\n" >> outputR "/" seqname[2] ".fa"}}' $genome

packageName="BSgenome."$organism"."$source".1"

echo "Package: "$packageName > seedFile.txt
echo "Title: Full genome sequences for "$organism >>  seedFile.txt
echo "Description: Full genome sequences for "$organism " from " $source >>  seedFile.txt
echo "Version: 1.0" >>  seedFile.txt
echo "organism: "$organism >>  seedFile.txt
echo "common_name: " $organism >>  seedFile.txt
echo "provider: " $source >>  seedFile.txt
echo "provider_version: 1" >>  seedFile.txt
echo "release_date: Apr. 2011" >>  seedFile.txt
echo "organism_biocview: "$organism >>  seedFile.txt
echo "release_name: Genome Sequencing" >>  seedFile.txt
echo "BSgenomeObjname: "$organism >>  seedFile.txt
echo "seqnames: "$namesSeq >>  seedFile.txt
echo "seqs_srcdir: "$outputRepository >>  seedFile.txt

Rscript /CirComPara/myScript/buildPackage.R

R CMD build $packageName
R CMD check $packageName"_1.0.tar.gz"
R CMD INSTALL $packageName"_1.0.tar.gz"






