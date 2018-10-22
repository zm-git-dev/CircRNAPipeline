library("dplyr");




  
  jsonConfig = snakemake@config
  
  dataAnalyseCirc = data.frame()
  dataAnalyseProcess = data.frame()
  dataHTSeqCount = data.frame()
  
  for(SampleIndex in 1:length(jsonConfig$SAMPLE)){
    SampleName = jsonConfig$SAMPLE[SampleIndex]
    
    Sample = jsonConfig[[SampleName]]
    pathCirc = paste0(Sample$REPLICAT, "/", Sample$CONDITION_NAME, "/Global_results/", SampleName, Sample$CONDITION_NAME, "_Annotate.csv")
    pathProcess = paste0(Sample$REPLICAT, "/", Sample$CONDITION_NAME, "/Global_results/", SampleName, Sample$CONDITION_NAME,"_", as.character(jsonConfig$USED_READS),"_ProcessingInfos.csv")
    pathHTSeqCount = paste0(Sample$REPLICAT, "/", Sample$CONDITION_NAME, "/Global_results/", SampleName, Sample$CONDITION_NAME,"_", as.character(jsonConfig$USED_READS_GENES), "_htSeqCount.txt")
    
    # On met tout sur la forme d'un tableau
    CircPrediction = read.table(file = pathCirc, header = TRUE, sep="\t", comment.char = "")
    CircPrediction$Sample = SampleName
    CircPrediction$Condition = Sample$CONDITION_NAME
    CircPrediction$Replicat = Sample$REPLICAT
    
    # ProcessInfo : On additionne les tables si on a utilisé à la fois les merge et les unmerged
    ProcessInfo = NULL
    if(length(pathProcess) == 2){
          ProcessInfo = read.table(file = pathProcess[1], header = TRUE, sep="\t") + read.table(file = pathProcess[2], header = TRUE, sep="\t")
    }else{
          ProcessInfo = read.table(file = pathProcess, header = TRUE, sep="\t")
    }
    ProcessInfo$Sample = SampleName
    ProcessInfo$Condition = Sample$CONDITION_NAME
    ProcessInfo$Replicat = Sample$REPLICAT
    
    # HTSeqCount : On additionne les matrices si on a utilisé à la fois les merge et les unmerged
    HTSeqCount = NULL
    if(length(pathHTSeqCount) == 2){
          HTSeqCount = read.table(file = pathHTSeqCount[1], sep = "\t", row.names = 1) + read.table(file = pathHTSeqCount[2], sep = "\t", row.names = 1)
    }else{
          HTSeqCount = read.table(file = pathHTSeqCount, sep = "\t", row.names = 1)
    }
    colnames(HTSeqCount) = SampleName
    
    dataAnalyseCirc = rbind(dataAnalyseCirc, CircPrediction)
    dataAnalyseProcess = rbind(dataAnalyseProcess, ProcessInfo)
    
   if(length(dataHTSeqCount) == 0){
	dataHTSeqCount = data.frame(GeneID = rownames(HTSeqCount), HTSeqCount)
    }else{
	dataHTSeqCount = cbind(dataHTSeqCount, HTSeqCount)
    }
    
    
  }
  saveRDS(object = dataAnalyseCirc, file = "./dataAnalyseCirc.rds")
  saveRDS(object = dataAnalyseProcess, file = "./dataAnalyseProcess.rds")
  saveRDS(object = dataHTSeqCount, file = "./dataHTSEQcount.rds")

