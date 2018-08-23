library("optparse");
library("magrittr");
library("rtracklayer");
library("biomaRt");
library("GenomicRanges")

##COMMAND EXAMPLE : Rscript RannotateFromGTF.R --mergeTable="./Janvier2014/Senescence/SenescenceFinalMergeTest.csv" --condition="Senescence" --path="./Janvier2014/" --gtf="./newGTF.gtf" --biomaRtDataSet="hsapiens_gene_ensembl"

##ATTENTION : Il FAUT rajouter une opt pour indiquer le path de output et l'intégrer dans le write.table sinon les règles de snakemake ne fonctioneront pas

getArguements = function(){
  option_list = list(
    make_option(c("--mergeTable"), type="character", default=NULL, 
                help=" Final with merged results from differents methods", metavar="character"),
    make_option(c("--condition"), type = "character", default = NULL,
                help="condition name", metavar="character"),
    make_option(c("--gtf"), type = "character", default = NULL,
                help="path of the GTF file", metavar="character"),
    make_option(c("--path"), type = "character", default = NULL,
                help="path of the ouput file", metavar="character")
  ); 
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  return(opt);
}



opt = getArguements();
condition = as.character(opt$condition);
path = as.character(opt$path);
mergeTable = as.character(opt$mergeTable);
gtf = as.character(opt$gtf);


###################################################################################################
################################### START FONCTIONS  ##############################################
###################################################################################################


TableOfClass = list("1" = "Exonic", "2" = "Intronic", "11" = "AntisensExonic", "12" = "AntisensIntronic", "0" = "Intergenic")

newAnnotateGR = function(GR, test, strand, resultat, AllGenesM){
  # On recherche tout les gènes mappant peut importe le brin afin de connaitre TOUT les transcripts et tout les gènes, a la fois pour le Start et pour le End tel que l'annotation puisse se faire par gènes
  Genes = genes[genes$gene_id %in% AllGenesM, ]
  # Si il y a au moins un gènes qui mappe sur le start et/ou sur le end alors on annote, sinon il s'agit d'un intergenique
  if(length(Genes) > 0){
  # On récupère tout les transcripts
  Tx = transcripts[subjectHits(findOverlaps(GR, transcripts)),]
  
  resultat[[paste0("TranscriptIds", test)]] = Tx$transcript_id
  V_ElementLocus = V_Class = vector("character", length(Genes))
  V_ExactMatchTxIds = V_ELementTxGlobal = list()
  
      # Pour chaque gène : DANS L'ORDRE
      for(i in (1:length(Genes))){
          gene = Genes[i,]
          PreviousV_Element = NULL
          # On cherche les élements génomique en Strand et Un_Strand
          if(strand == as.character(strand(gene))){
            PreviousV_Element = computeAnnotation(GR, gene$gene_id, exonsByGenes, intronsInGenes, "@E", "@I", as.character(strand(gene)), strand, 0)
          }else{
            PreviousV_Element = computeAnnotation(GR, gene$gene_id, exonsByGenes, intronsInGenes, "@AE", "@AI", as.character(strand(gene)), ifelse(strand == "+", "-", "+"), 10)
          }
          # En fait pour le gène i au strat ou au end  si il y a un élément génomique il n' en a qu'un donc on regarde si la concat est > 0 si non c'est donc que rien n'as mappé on met NA
          V_ElementLocus[i] = ifelse(length(PreviousV_Element$Element) > 0, yes = PreviousV_Element$Element, no = "NA")
          V_Class[i] = PreviousV_Element$Class
          
          TxForGene_i = Tx[Tx$gene_id == gene$gene_id,] 
          if(length(TxForGene_i) > 0){
              V_ELementTx = rep("", length(TxForGene_i))
              for(j in (1:length(TxForGene_i))){
                transcript = TxForGene_i[j, ]
                PreviousV_Element = NULL
                if(strand == as.character(strand(gene))){
                  PreviousV_Element = computeAnnotation(GR, transcript$transcript_id, exonsByTranscripts, intronsInTranscripts, "@E", "@I", as.character(strand(transcript)), strand, 0)$Element
                }else{
                  PreviousV_Element = computeAnnotation(GR, transcript$transcript_id, exonsByTranscripts, intronsInTranscripts, "@AE", "@AI", as.character(strand(transcript)), ifelse(strand == "+", "-", "+"), 10)$Element
                }
                V_ELementTx[j] = ifelse(length(PreviousV_Element) > 0, yes = PreviousV_Element, no = "NA")
              }
              V_ELementTxGlobal[[i]] = paste0(V_ELementTx, collapse = ",")
          }else{
            V_ELementTxGlobal[[i]] = "NA"
          }
          
          # On fait le test de l'exact match à l'échelle du gène surlequel on travaille
          TxKeepedExactMatchBool = testIfExactMatch(TxForGene_i, GR, test)
          V_ExactMatchTxIds[[i]] = TxForGene_i[TxKeepedExactMatchBool, ]$transcript_id
      }
      V_Class = as.vector(unlist(TableOfClass[as.character(V_Class)]))
      
  resultat[[paste0("ElementLocus",test)]] = V_ElementLocus
  resultat[[paste0("ElementTranscripts",test)]] = unlist(V_ELementTxGlobal)
  resultat[[paste0("Class",test)]] = V_Class
  resultat[[paste0("ExactMatchTxIds",test)]] = V_ExactMatchTxIds
  
  }else{
    return(annotateIntergenic(resultat, test))
  }
  
  return(resultat)
}


# idToBrowse: Les ids de gènes ou de Tx
# ExonsByX: exonsBygenes ou exonsBytx
# IntronsByX: pareil
# sepForExons et sepForIntrons: les @...
#: GenomicElementTable Genes ou Transcripts
computeAnnotation = function(GR, idToBrowse, ExonsByX, IntronsByX, sepForExons, sepForIntrons, StrandGenomicElement, strand, factor){
  # On utilise un code pour la classe avec intergenic a 0 pour faciliter l'interprétation (0 + 0) ou 0+x donne toujorus x vue que l'on appelle deux fois la méthode et que elle ne mappera au max que a un endroit
  V_New_Element = NULL
  V_Class = 0
  strand(GR) = strand
    # On récupère les éléments G associés
    Exons = sort(ExonsByX[[idToBrowse]] %>% reduce(), decreasing = (StrandGenomicElement == "-"))
    introns = IntronsByX[[idToBrowse]]
    # Si il il est associé à des exons alors on envoie
    TestExon = subjectHits(findOverlaps(GR, Exons))
    if(length(TestExon) > 0){
      V_New_Element = paste(idToBrowse, TestExon, collapse = ",", sep = sepForExons )
      V_Class = V_Class + 1 + factor
    }else if(!is.null(introns)){
      # Car un exon peut être constitué d'un seul exons et ne pas avoir d'introns
      TestIntron = subjectHits(findOverlaps(GR,introns))
      if(length(TestIntron) > 0){
        V_New_Element = paste(idToBrowse, introns[TestIntron,]$intron_id , collapse = ",", sep = sepForIntrons )
        V_Class = V_Class + 2 + factor
      }
    }
  return(list(Element = V_New_Element, Class = V_Class))
}


# Annota est la méthode qui va permettre d'obtenir une annotation en start et en End (on ne peut pas vraiment savoir ce qu'il y a au milieu du circRNA) ensuite on merge les résultat des annot start et end
annotate = function(dataLine, annotateData){
  start = dataLine$Start
  end = dataLine$End
  chr = as.character(dataLine$Chromosome)
  strand = as.character(dataLine$Strand)
  
  resultat = list();
  # On crée des GR initialement sans brin (pour les cas antisens)
  GRStartWithoutStrand = GRanges(seqnames = chr, ranges = IRanges(start+1, start+2+confidenceWindow),strand = "*")
  GREndWithoutStrand = GRanges(seqnames = chr, ranges = IRanges(end-1-confidenceWindow, end),strand = "*");
  #On récupère les Genes en Start et les Genes et End
  GenesMatchInStart = genes[subjectHits(findOverlaps(GRStartWithoutStrand, genes)),]$gene_id
  GenesMatchInEnd = genes[subjectHits(findOverlaps(GREndWithoutStrand, genes)),]$gene_id
  AllGenesM = base::union(GenesMatchInStart, GenesMatchInEnd)
  resultat = newAnnotateGR(GRStartWithoutStrand, "Start", strand, resultat, AllGenesM);
  resultat = newAnnotateGR(GREndWithoutStrand, "End", strand, resultat, AllGenesM);
  resultat[is.na(resultat)] = ""
  Genes = genes[genes$gene_id %in% AllGenesM, ]
  
  nbGenes = GenesMatchInStartNoAntisens = GenesMatchInEndNoAntisens = SelectedGenesNoAntinsens = NULL
  # Pour détermines les chimères, on ne considère pas les transcripts antisens car ils pourrait créer de fausses chimères
  if(length(AllGenesM) > 0){
    nbGenes = length(AllGenesM)
    SelectedGenesNoAntinsens = AllGenesM[((resultat[["ClassStart"]] != "AntisensExonic" & resultat[["ClassStart"]] != "AntisensIntronic") & (resultat[["ClassEnd"]] != "AntisensExonic" & resultat[["ClassEnd"]] != "AntisensIntronic"))]
    GenesMatchInStartNoAntisens = GenesMatchInStart[GenesMatchInStart %in% SelectedGenesNoAntinsens]
    GenesMatchInEndNoAntisens = GenesMatchInEnd[GenesMatchInEnd %in% GenesMatchInEnd]
  }else{
    nbGenes = 1
  }
  
  
  # Si il a été annoté a l'échelle du gène comme étant intergenic mais que il y a au moins un gène qui mappe alors il est pas vraiment intergenic, juste notIncluded dans ce gène
  for(i in 1:nbGenes){
    if(resultat[["ClassStart"]][i] == "Intergenic" && length(GenesMatchInStart) > 0){
      resultat[["ClassStart"]][i] = "notIncluded"
    }else if(resultat[["ClassEnd"]][i] == "Intergenic" && length(GenesMatchInEnd) > 0){
      resultat[["ClassEnd"]][i] = "notIncluded"
    }
  }
  
  #On initiale tout a false, pour les intergenic et out c'est plus pratique
  isChimeraList = rep(FALSE, nbGenes)
  # Pour un gènes qui n'est pas un antisens: Si il y a au moins 1 gènes (non antisens) en start et en end ET que l'intersection de ces gènes est nulle (aucun communs) alors ces gènes (hotmis les antisens) sont des chimères
  if((length(GenesMatchInStartNoAntisens) > 0) & (length(GenesMatchInEndNoAntisens) > 0) ){
    isChimeraList[which(AllGenesM %in% SelectedGenesNoAntinsens)] = (length(intersect(GenesMatchInStartNoAntisens, GenesMatchInEndNoAntisens)) == 0)
  }
  
  newDataAnnotate = dataLine[rep(1, nbGenes),]
  newDataAnnotate$Class = resumeAnnotationFactor(resultat, "Class", "-")
  newDataAnnotate$TranscriptIds = resumeAllAnnotationFactor(resultat, "TranscriptIds", ",")
  newDataAnnotate$GeneId = paste0(Genes$gene_id, "")
  newDataAnnotate$GeneSymbol = paste0(Genes$gene_symbol, "")
  newDataAnnotate$GeneBioType = paste0(Genes$gene_biotype, "")
  newDataAnnotate$ElementLocusStart = resultat[["ElementLocusStart"]]
  newDataAnnotate$ElementLocusEnd = resultat[["ElementLocusEnd"]]
  newDataAnnotate$GlobalElementLocusStart = paste0(resultat[["ElementLocusStart"]], collapse = ",")
  newDataAnnotate$GlobalElementLocusEnd = paste0(resultat[["ElementLocusEnd"]], collapse = ",")
  newDataAnnotate$ElementTranscriptsStart = resultat[["ElementTranscriptsStart"]]
  newDataAnnotate$ElementTranscriptsEnd = resultat[["ElementTranscriptsEnd"]]
  newDataAnnotate$GlobalElementTranscriptsStart = paste0(resultat[["ElementTranscriptsStart"]][resultat[["ElementTranscriptsStart"]] != "NA"], collapse = ",")
  newDataAnnotate$GloablElementTranscriptsEnd = paste0(resultat[["ElementTranscriptsEnd"]][resultat[["ElementTranscriptsEnd"]] != "NA"], collapse = ",")
  
  newDataAnnotate$ExactMatch =  sapply(1:nbGenes, FUN = function(i){ ifelse( length(intersect(resultat[["ExactMatchTxIdsStart"]][[i]], resultat[["ExactMatchTxIdsEnd"]][[i]])) > 0, TRUE, FALSE) })
  newDataAnnotate$isChimera = isChimeraList
  newDataAnnotate$ExactMatchIds = sapply(1:nbGenes, FUN = function(i){concatRes(intersect(resultat[["ExactMatchTxIdsStart"]][[i]], resultat[["ExactMatchTxIdsEnd"]][[i]]), ",") })

  annotateData = rbind(annotateData, newDataAnnotate)
  
  return(annotateData);
}


#Ici on fait par élement pour les gènes distinct d'une même jinction: Class, gene_id, ...
resumeAnnotationFactor = function(resultat, factor, separateur){
  merging = c()
  for(i in 1:length(resultat[[paste0(factor,"Start")]])){
    merging = c(merging, concatRes(base::union(resultat[[paste0(factor,"Start")]][i], resultat[[paste0(factor,"End")]][i]), separateur))
  }
  return(merging)
  }

#Ici on fait le global peut importe donc on fait la totale, notamment pour Transcripts, ..
resumeAllAnnotationFactor = function(resultat, factor, separateur){
  merging = concatRes(base::union(resultat[[paste0(factor,"Start")]], resultat[[paste0(factor,"End")]]), separateur)
  return(merging)
}

concatRes = function(element, separateur){
  return(paste0(element, sep="", collapse = separateur));
}



testIfExactMatch = function(GROverlapTx, GRElement, test){
  #Si test = start alors on test le start de l'exon A sinonn le end de l'exon B dans le sens du génome
  positonQuery = if(test == "Start") start(GRElement) else end(GRElement)
  #On récupère tout les exons associés aux transcripts mappés en liste
  AllExonsOfSelectTranscripts = exonsByTranscripts[which(names(exonsByTranscripts) %in% names(GROverlapTx))]
  #suivant si Start ou End, on regarde lesquels sont égaux aux start ou end
  # ATTENTION vue que start en end peuvent faires des exacts match indépendant on supposer qu'il peut arriver que lors de l'intersection aucun ne soit commum, donc on peut imaginer ExactMatch True mais pas de ExactMatchIds...
  TxKeeped = if(test == "Start") lapply(AllExonsOfSelectTranscripts, FUN = function(x){ return( ifelse((positonQuery %in% start(x)), TRUE, FALSE) ) })  else lapply(AllExonsOfSelectTranscripts, FUN = function(x){ return( ifelse((positonQuery %in% end(x)), TRUE, FALSE) ) })
  return(unlist(TxKeeped, use.names = FALSE))
  
  }



#Annotation pour les position Intergéniques
annotateIntergenic = function(resultat, test){
  resultat[[paste0("Class",test)]] = "Intergenic";
  resultat[[paste0("TranscriptIds", test)]] = resultat[[paste0("GeneId",test)]] = resultat[[paste0("GeneSymbol",test)]] = resultat[[paste0("GeneBioType",test)]] = resultat[[paste0("ElementLocus",test)]] = resultat[[paste0("ElementTranscripts",test)]] = resultat[[paste0("ElementLocus", test)]] = resultat[[paste0("ElementTranscripts", test)]] = "";
  resultat[[paste0("ExactMatchTxIds", test)]] = c()
  return(resultat);
}




  ###################################################################################################
  #################################### END FONCTIONS  ###############################################
  ###################################################################################################


  ################################### GET_OBJECTS  ##################################################
  
  
  exons = readRDS("./Rsave/exons.rds");
  genes = readRDS("./Rsave/genes.rds");
  transcripts = readRDS("./Rsave/transcripts.rds");
  exonsByGenes = readRDS("./Rsave/exonsByGenes.rds");
  exonsByTranscripts = readRDS("./Rsave/exonsByTranscripts.rds");
  genesWithoutStrand = readRDS("./Rsave/genesWithoutStrand.rds");
  intronsInTranscripts = readRDS("./Rsave/intronsInTranscripts.rds");
  intronsInGenes = readRDS("./Rsave/intronsInGenes.rds");
  
  intronsInTranscripts = split(intronsInTranscripts, intronsInTranscripts$transcript_id); 
  intronsInGenes = split(intronsInGenes, intronsInGenes$gene_id);

  #####################################  ANNOTATION  ################################################


print("Begin Annotation of candidates ....") 

#On importe le GTF
input = import(gtf);

confidenceWindow = 0; 

chrList = as.character(levels(seqnames(input)))

data = read.table(mergeTable, header = TRUE, sep = "\t", comment.char = "");

completeAnnotation = function(data, chrList){
  
  annotatedData = data.frame()
  
  for(i in 1:nrow(data)){
    
    if(data[i,]$Chromosome %in% chrList){
      
      annotatedData = annotate(data[i,], annotatedData)
      
      }
  
    }
  return(annotatedData)
}

resultat = completeAnnotation(data, chrList)
write.table(resultat, file = paste0(path, condition, "_Annotate.csv"),row.names=FALSE, na="", sep = "\t", dec = ".", quote = FALSE)


