library("Biostrings");
library("GenomicFeatures");
library("rtracklayer");
library("dplyr")
library("lazyeval")

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

condition = as.character(snakemake@params$condition)
print(condition)
path = as.character(snakemake@params$path)
print(path)
packageNameGenome = as.character(snakemake@params$packageName)
print(packageNameGenome)

pathCircExplorer2 = as.character(snakemake@input$CIRCEXPLORER2)
print(pathCircExplorer2)
pathCircRNA_finder = as.character(snakemake@input$CIRCRNAFINDER)
print(pathCircRNA_finder)
pathFind_circ = as.character(snakemake@input$FINDCIRC)
print(pathFind_circ)

# Load Pipemade package for genome
library(packageNameGenome, character.only = TRUE)
#Variable globale
genome <- getBSgenome(packageNameGenome);


getProcessedFile = function(ListOfFiles, ToolName){
  ProcessedTable = NULL;
  switch(ToolName,
         CircExplorer2={
           ProcessedTable = getPre_processedFile(ListOfFiles, "Junctions_reads_CircExplorer2", "Reads_CircExplorer2", "PreProcess_CricExplorer2")
           ProcessedTable = getLength(ProcessedTable);
           ProcessedTable = catchOnlySignal(ProcessedTable);
         },
         circRNAFinder={
           ProcessedTable = getPre_processedFile(ListOfFiles, "Junctions_reads_circRNA_finder", "Reads_circRNAFinder", "PreProcess_CircRNAFinder")
           ProcessedTable = getLength(ProcessedTable);
           ProcessedTable = catchOnlySignal(ProcessedTable);
         },
         findCirc={
           ProcessedTable = getPre_processedFile(ListOfFiles, "Junctions_reads_FindCirc", "Id_Read_FindCirc", "PreProcess_Find_circ")
           ProcessedTable = getLength(ProcessedTable);
         }
  )
}

getPre_processedFile = function(ListOfFiles, junctionFactor, readFactor, PreProcessFunction){
  # Si on a fait merge et unmerge
  Pre_processedTable = NULL;
  if(length(ListOfFiles) == 2){
      Pre_processedTable = mergePrediction(do.call(PreProcessFunction, list(ListOfFiles[1])), do.call(PreProcessFunction, list(ListOfFiles[2])), junctionFactor, readFactor )
  }else{
      Pre_processedTable = do.call(PreProcessFunction, list(ListOfFiles[1]))
  }
  return(as.data.frame(Pre_processedTable))
}

mergePrediction = function(dataFromMerged, dataFromUnmerged, junctionFactor, readFactor){
  
  mutate_call_junction = lazyeval::interp(~ sum(as.numeric(junctions)), junctions = as.name(junctionFactor))
  mutate_call_reads = lazyeval::interp(~ paste(reads, collapse = ""), reads = as.name(readFactor))
  
  mergedData = rbind(dataFromMerged, dataFromUnmerged) %>% group_by(Chromosome, Start, End, Strand) %>% mutate_(.dots = setNames(list(mutate_call_junction), junctionFactor)) %>% mutate_(.dots = setNames(list(mutate_call_reads), readFactor)) %>% distinct()
  return(mergedData)
}

PreProcess_CricExplorer2 = function(file_CircExplorer2){
  
      if(!is.null(file_CircExplorer2)){
          dataCircExplorer2 = read.table(file_CircExplorer2, header = TRUE, sep = "\t", comment.char = "");
          dataCircExplorer2 = subset.data.frame(dataCircExplorer2, select = c(Chromosome, Start_of_fusion_junction, End_of_fusion_junction, Fusion_id.Junction_reads, Strand, Reads))
          cc<- strsplit(as.character(dataCircExplorer2$Fusion_id.Junction_reads),"/")
          dataCircExplorer2$Fusion_id.Junction_reads = unlist(cc)[2*(1:length(dataCircExplorer2$Fusion_id.Junction_reads)) ]
          colnames(dataCircExplorer2) = c("Chromosome", "Start", "End", "Junctions_reads_CircExplorer2", "Strand", "Reads_CircExplorer2");
          return(dataCircExplorer2)
      }else{
          return(NULL)
      }
}
 
PreProcess_CircRNAFinder = function(file_CircRNAFinder){
        if(!is.null(file_CircRNAFinder)){
            dataCircRNAFinder = read.table(file_CircRNAFinder, header = TRUE, sep = "\t", comment.char = "");
            #On inverse les brin car dans le script ils les retournent à l'envers donc bon ....
            dataCircRNAFinder$Strand = apply(as.matrix(dataCircRNAFinder$Strand), 1, function(x){if(as.character(x) == "+") "-" else "+"})
            dataCircRNAFinder = subset.data.frame(dataCircRNAFinder, select = c(Chromosome, Start, End, Number_of_junction_reads, Strand, Reads));
            colnames(dataCircRNAFinder) = c("Chromosome", "Start", "End", "Junctions_reads_circRNA_finder", "Strand","Reads_circRNAFinder");
            return(dataCircRNAFinder)
        }else{
          return(NULL)
        }
}

PreProcess_Find_circ = function(file_findCirc){
        dataFindCirc = read.table(file_findCirc, header = TRUE, sep = "\t", comment.char = "");
        dataFindCirc = subset.data.frame(dataFindCirc, select = c(Chromosome, left_splice_site, right_splice_site, number_of_junction_reads, strand, signal, name));
        colnames(dataFindCirc) = c("Chromosome", "Start", "End", "Junctions_reads_FindCirc", "Strand", "Signal", "Id_Read_FindCirc");
        return(dataFindCirc)
}

dataCircExplorer2 = getProcessedFile(pathCircExplorer2, "CircExplorer2")
dataCircRNAFinder = getProcessedFile(pathCircRNA_finder, "circRNAFinder")
dataFindCirc = getProcessedFile(pathFind_circ, "findCirc")


firstMerge = merge.data.frame(x = dataFindCirc, y =  dataCircRNAFinder, by.x = c("Chromosome", "Start", "End", "Strand", "Length", "Signal"), by.y = c("Chromosome", "Start", "End", "Strand", "Length", "Signal"), all = TRUE);

FinalMerge = merge.data.frame(x = firstMerge, y = dataCircExplorer2, by.x = c("Chromosome", "Start", "End", "Length", "Strand", "Signal"), by.y = c("Chromosome", "Start", "End", "Length", "Strand", "Signal"), all = TRUE);

FinalMerge$supportedMethod = apply(subset.data.frame(FinalMerge, select = c(Junctions_reads_CircExplorer2, Junctions_reads_circRNA_finder, Junctions_reads_FindCirc)),1,function(x){sum(!is.na(x))})


write.table(FinalMerge, file = paste0(path, condition, "Merge.csv"),row.names=FALSE, na="", sep = "\t", dec = ".", quote = FALSE)





