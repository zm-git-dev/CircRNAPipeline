library("optparse");
# library("BSgenome.Hsapiens.UCSC.hg38"); 
library("Biostrings");
library("GenomicFeatures");
library("rtracklayer");


#COMANNDE : Rscript mergeScript.R --circExplorer2="./HMKMGBGXY_Sen_back_spliced_junction_withHeaders.bed" --circRNAFinder="HMKMGBGXY_Sen_filteredJunctions_withHeaders.bed" --findCirc="HMKMGBGXY_Sen_filtered_splice_circ_withHeaders.bed" --condition="Senescence" --outputPath="./Octobre2016/"

getArguements = function(){
option_list = list(
  make_option(c("--circExplorer2"), type="character", default=NULL, 
              help=" backspliced file name for circExplorer2", metavar="character"),
  make_option(c("--circRNAFinder"), type="character", default=NULL, 
              help="circRNA candidates from circRNAFinder", metavar="character"),
  make_option(c("--findCirc"), type = "character", default = NULL,
        help="circRNA candidates from findCirc", metavar="character"),
  make_option(c("--condition"), type = "character", default = NULL,
              help="circRNA candidates from findCirc", metavar="character"),
  make_option(c("--outputPath"), type = "character", default = NULL,
              help="Path for output file", metavar="character"),
  make_option(c("--packageName"), type = "character", default = NULL,
              help="packageName for genome", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
return(opt);
}


getLength = function(data){
  data$Length = data$End - data$Start;
  return(data);
}

catchOnlySignal = function(data){
  
  start = getSeq(genome, as.vector(data$Chromosome), start = as.vector(data$Start) -1, end = as.vector(data$Start));
  end = getSeq(genome, as.vector(data$Chromosome), start = as.vector(data$End) +1, end = as.vector(data$End) +2);
  site = paste0(end, start);
  finalsite = sapply(1:nrow(data),function(x){ifelse(data[x, ]$Strand == "+", as.character(site[x]), as.character(reverseComplement(DNAString(site[x]))))});
  data$Signal = finalsite;
  return(data);
  
}


########
# MAIN #
########


opt = getArguements();

condition = as.character(opt$condition);

# condition = "Senescence"
path = as.character(opt$outputPath);
packageNameGenome = as.character(opt$packageName);
# Load Pipemade package for genome
library(packageNameGenome, character.only = TRUE)

#Variable globale
genome <- getBSgenome(packageNameGenome);

dataCircExplorer2 = read.table(opt$circExplorer2, header = TRUE, sep = "\t", comment.char = "");
# dataCircExplorer2 = read.table("./Janvier2014/Senescence/SEN_L1_back_spliced_junction_withHeaders.bed", header = TRUE, sep = "\t");


#Select column from CircExplorer2 and Treatment
dataCircExplorer2 = subset.data.frame(dataCircExplorer2, select = c(Chromosome, Start_of_fusion_junction, End_of_fusion_junction, Fusion_id.Junction_reads, Strand, Reads))
cc<- strsplit(as.character(dataCircExplorer2$Fusion_id.Junction_reads),"/")
dataCircExplorer2$Fusion_id.Junction_reads = unlist(cc)[2*(1:length(dataCircExplorer2$Fusion_id.Junction_reads)) ]
colnames(dataCircExplorer2) = c("Chromosome", "Start", "End", "Junctions_reads_CircExplorer2", "Strand", "Reads_CircExplorer2");
dataCircExplorer2 = getLength(dataCircExplorer2);
dataCircExplorer2 = catchOnlySignal(dataCircExplorer2);



dataCircRNAFinder = read.table(opt$circRNAFinder, header = TRUE, sep = "\t", comment.char = "");
# dataCircRNAFinder = read.table("./Janvier2014/Senescence/SEN_L1_filteredJunctions_withHeaders.bed", header = TRUE, sep = "\t");

#On inverse les brin car dans le script ils les retournent à l'envers donc bon ....
#Select column from CircExplorer2 and Treatment
dataCircRNAFinder$Strand = apply(as.matrix(dataCircRNAFinder$Strand), 1, function(x){if(as.character(x) == "+") "-" else "+"})
dataCircRNAFinder =subset.data.frame(dataCircRNAFinder, select = c(Chromosome, Start, End, Number_of_junction_reads, Strand, Reads));
colnames(dataCircRNAFinder) = c("Chromosome", "Start", "End", "Junctions_reads_circRNA_finder", "Strand","Reads_circRNAFinder");
dataCircRNAFinder = getLength(dataCircRNAFinder);
dataCircRNAFinder = catchOnlySignal(dataCircRNAFinder);



dataFindCirc = read.table(opt$findCirc, header = TRUE, sep = "\t", comment.char = "");

#Select column from FindCirc and Treatment
dataFindCirc =subset.data.frame(dataFindCirc, select = c(Chromosome, left_splice_site, right_splice_site, number_of_junction_reads, strand, signal, name));
colnames(dataFindCirc) = c("Chromosome", "Start", "End", "Junctions_reads_FindCirc", "Strand", "Signal", "Id_Read_FindCirc");
dataFindCirc = getLength(dataFindCirc);


firstMerge = merge.data.frame(x = dataFindCirc, y =  dataCircRNAFinder, by.x = c("Chromosome", "Start", "End", "Strand", "Length", "Signal"), by.y = c("Chromosome", "Start", "End", "Strand", "Length", "Signal"), all = TRUE);

FinalMerge = merge.data.frame(x = firstMerge, y = dataCircExplorer2, by.x = c("Chromosome", "Start", "End", "Length", "Strand", "Signal"), by.y = c("Chromosome", "Start", "End", "Length", "Strand", "Signal"), all = TRUE);

FinalMerge$supportedMethod = apply(subset.data.frame(FinalMerge, select = c(Junctions_reads_CircExplorer2, Junctions_reads_circRNA_finder, Junctions_reads_FindCirc)),1,function(x){sum(!is.na(x))})


#write.csv(FinalMerge, file = "FinalMergeTest.csv",row.names=FALSE, na="")
write.table(FinalMerge, file = paste0(path, condition, "Merge.csv"),row.names=FALSE, na="", sep = "\t", dec = ".", quote = FALSE)





