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
    make_option(c("--path"), type = "character", default = NULL,
                help="path of the ouput file", metavar="character"),
    make_option(c("--gtf"), type = "character", default = NULL,
                help="path of the GTF file", metavar="character"),
    make_option(c("--biomaRtDataSet"), type = "character", default = NULL,
                help="Name of the biomaRt dataset to use for annotation", metavar="character")
  ); 
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  return(opt);
}



opt = getArguements();
condition = as.character(opt$condition);
path = as.character(opt$path);
gtf = as.character(opt$gtf);
mergeTable = as.character(opt$mergeTable);
# biomaRtDataSet est directement utilisé dans le if pour savoir s'il a été envoyé sinon il n'y aura pas de geneSymbol et GeneBioType

#### Avec RtrackLayer

# gtf = "./newGTF.gtf";
# path = "./Janvier2014/";
# condition ="Senescence";
# biomaRtDataset = "hsapiens_gene_ensembl";
# mergeTable = "./Janvier2014/Senescence/SenescenceFinalMergeTest.csv";
 

###################################################################################################
################################### START FONCTIONS  ##############################################
###################################################################################################



# Annota est la méthode qui va permettre d'obtenir une annotation en start et en End (on ne peut pas vraiment savoir ce qu'il y a au milieu du circRNA) ensuite on merge les résultat des annot start et end
annotate = function(start, end, chr, strand){
  resultat = list();
  # On crée des GR initialement sans brin (pour les cas antisens)
  GRStartWithoutStrand = GRanges(seqnames = chr, ranges = IRanges(start+1, start+2+confidenceWindow),strand = "*")
  GREndWithoutStrand = GRanges(seqnames = chr, ranges = IRanges(end-1-confidenceWindow, end),strand = "*");
  resultat = annotateGR(GRStartWithoutStrand, "Start", strand, resultat);
  resultat = annotateGR(GREndWithoutStrand, "End", strand, resultat);
  resultat[is.na(resultat)] = ""
  resultat = mergeResultat(resultat);
  return(resultat);
}

concatRes = function(element){
  return(paste0(element, sep="", collapse = ','));
}


mergeResultat = function(resultat){
  resultat = getCommonAnnotationForList(resultat, "TranscriptIds");
  resultat = getCommonAnnotation(resultat, "Class");
  resultat = getCommonAnnotation(resultat, "GeneId");
  resultat = getCommonAnnotation(resultat, "GeneSymbol");
  resultat = getCommonAnnotation(resultat, "GeneBioType");
  resultat = getCommonAnnotationForList(resultat, "ExactMatchTxIds");
  resultat[["ExactMatch"]] = resultat[["ExactMatchStart"]] & resultat[["ExactMatchEnd"]]
  return(resultat);
}

# On recherche les annotation communes entre start et end pour éviter la redondance
getCommonAnnotationForList = function(resultat, factor){
  if(!setequal(resultat[[paste0(factor,"Start")]], resultat[[paste0(factor,"End")]])){
    if(resultat[["GeneIdStart"]] == resultat[["GeneIdEnd"]]){
      resultat[[factor]] = concatRes(base::intersect(resultat[[paste0(factor,"Start")]], resultat[[paste0(factor,"End")]])) #Si il s'agit de même gène alors on fait l'intersection
      # resultat[[factor]] = concatRes(unique(append(resultat[[paste0(factor,"Start")]], resultat[[paste0(factor,"End")]])))
    }else{
      resultat[[factor]] = paste0("Start:", concatRes(resultat[[paste0(factor,"Start")]]),"||End:", concatRes(resultat[[paste0(factor,"End")]]));
    } 
  }else{
    resultat[[factor]] = concatRes(resultat[[paste0(factor,"Start")]]);
  }
  return(resultat);
}

# On recherche les annotation communes entre start et end pour éviter la redondance
getCommonAnnotation = function(resultat, factor){
  if(resultat[[paste0(factor,"Start")]] != resultat[[paste0(factor,"End")]]){
    resultat[[factor]] = paste0("Start:", resultat[[paste0(factor,"Start")]],"||End:", resultat[[paste0(factor,"End")]]);
  }else{
    resultat[[factor]] = resultat[[paste0(factor,"Start")]];
  }
  return(resultat);
}

# Permet de réaliser une suite de Tests pour annoter Si on mappe un Gène sans brin alors il est soit "exonic" ou "intronic" ou "antisens" sinon c'est "intergenic"
annotateGR = function(GR, test, strand, resultat){
  GeneOverloap = findOverlaps(GR, AllTranscriptWithoutStrandWithIds);
  if(length(GeneOverloap) > 0){
    strand(GR) = strand;
    TElementTest = findOverlaps(GR, exonsByTranscripts);
    TElement = transcripts[which(transcripts$transcript_id %in% names(exonsByTranscripts[subjectHits(TElementTest), ])), ];
    IntronsTest = findOverlaps(GR, introns);
    if(length(TElement) != 0){
      resultat = annotateTranscript(TElement, GR, strand, resultat, test);
      return(resultat);
    }else if(length(IntronsTest) > 0){
      resultat = annotateIntrons(GR, IntronsTest, strand, resultat, test);
      return(resultat);
    }else{
      resultat =  annotateAntinsens(GR, strand, resultat, test);
      return(resultat)
    }
    
  }else{
    resultat = annotateIntergenic(resultat, test);
    return(resultat);
  }
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
  resultat[[paste0("TranscriptIds", test)]] = resultat[[paste0("GeneId",test)]] = resultat[[paste0("GeneSymbol",test)]] = resultat[[paste0("GeneBioType",test)]] = resultat[[paste0("Rank",test)]] = resultat[[paste0("Element",test)]] = resultat[[paste0("ExactMatchTxIds", test)]] = "";
  resultat[[paste0("ExactMatch", test)]] = FALSE
  return(resultat);
}

# Annotation pour les antisens
annotateAntinsens = function(GRGene, strand, resultat, test){
        resultat = annotateGR(GRGene, test, ifelse(strand == "+", "-", "+"), resultat)
        resultat[[paste0("Class",test)]] = "Antinsens";
        return(resultat)
}

#Annotation pour les introns
# ATTENTION a cause des pseudogène (appartenant a un même gène mais situé a des positions très différente il perturbe et donne de fausse reconstruction de Locus) donc on va tout de même vérifier pour les introns que l'on match bien entre le start et le end d'un transcrit
annotateIntrons = function(GR, IntronsTest, strand, resultat, test){
  #On révupère tout les introns matchant
  myIntrons = introns[subjectHits(IntronsTest), ];
    geneID = myIntrons$gene_id[1];
    resultat[[paste0("GeneId",test)]] = geneID
    resultat[[paste0("TranscriptIds", test)]] = transcripts[subjectHits(findOverlaps(GR, transcripts)), ]$transcript_id;
    resultat[[paste0("Class",test)]] = "Intronic";
    resultat[[paste0("GeneSymbol",test)]] = genes[genes$gene_id == geneID]$gene_symbol
    resultat[[paste0("GeneBioType",test)]] = genes[genes$gene_id == geneID]$gene_biotype
    resultat[[paste0("Rank",test)]] = myIntrons$intron_id;
    resultat[[paste0("Element",test)]] =  paste0("Intron", myIntrons$intron_id);
    resultat[[paste0("ExactMatch", test)]] = FALSE
    resultat[[paste0("ExactMatchTxIds", test)]] = ""
    return(resultat);
}

# Annotation pour les transcrit sur des match exoniques
annotateTranscript = function(TElement, GRElement, strand, resultat, test){
  #Doc Cf annotateAntisens
  TxKeepedExactMatchBool = testIfExactMatch(TElement, GRElement, test)
  TxExactMatch = TElement[TxKeepedExactMatchBool, ]
  resultat[[paste0("ExactMatch", test)]] = ifelse(sum(TxKeepedExactMatchBool) != 0, TRUE, FALSE)
  resultat[[paste0("ExactMatchTxIds", test)]] = TxExactMatch$transcript_id
  
  resultat[[paste0("TranscriptIds", test)]] = TElement$transcript_id
  geneID = as.character(TElement[1,]$gene_id[1]);
  resultat[[paste0("GeneId",test)]] = geneID;
  resultat[[paste0("GeneSymbol",test)]] = genes[genes$gene_id == geneID]$gene_symbol
  resultat[[paste0("GeneBioType",test)]] = genes[genes$gene_id == geneID]$gene_biotype
  resultat[[paste0("Class",test)]] = "Exonic";
  AllExons = reduce(exonsByGenes[[geneID]]);
  if(strand == "-"){
    AllExons = sort(AllExons, decreasing = TRUE);
  }else{
    AllExons = sort(AllExons); 
  }
  resultat[[paste0("Rank",test)]] = subjectHits(findOverlaps(GRElement, AllExons));
  resultat[[paste0("Element",test)]] = paste0("Exon", resultat[[paste0("Rank",test)]]);
  return(resultat);
}


# Crée un objet GRanges représentant les introns à paritr des coordonnées des exons de chaque Locus
getIntrons  = function(exonsByGenes, genes){
  GRIntrons = GRanges();
  for(geneName in names(exonsByGenes)){
    exons <- exonsByGenes[[geneName]] %>% reduce() %>% sort();
    if(length(exons) > 1){
      currentGRgene = GRanges(seqnames=seqnames(exons)[1], ranges=IRanges(start=min(start(exons)), end=max(end(exons))), strand=strand(exons)[1])  
      db = disjoin(c(exons, currentGRgene));
      ints = db[countOverlaps(db, exons) == 0]
      #Add an ID
      if(length(ints) > 0){
        if(as.character(strand(currentGRgene)) == "-") {
          ints$intron_id = c(length(ints):1)
        } else {
          ints$intron_id = c(1:length(ints))
        }
        
        ints$gene_id = geneName
        ints$symbol = genes[genes$gene_id == geneName,]$gene_symbol;
        GRIntrons = append(GRIntrons, ints);
      }  
    }
  }
  return(GRIntrons);
}


separateTxFromStAndChr = function(genes){
  
  for(geneId in names(genes)){
    currentGene = genes[[geneId]]
    if( length(unique(strand(currentGene))) > 1  || length(unique(seqnames(currentGene))) > 1){
      
      geneFilterOnChr = GRangesList(split(currentGene, seqnames(currentGene)))
      geneFilterOnChr = lapply(geneFilterOnChr, FUN = function(x){ if(length(x) > 0){ return(x) } } )
      geneFilterOnChr = unlist(geneFilterOnChr)
      geneFilterOnChr = lapply(X = geneFilterOnChr, function(x){ 
        geneFilterOnStrand = GRangesList(split(x, strand(x)))
        geneFilterOnStrand = lapply(geneFilterOnStrand, FUN = function(x){ if(length(x) > 0){ return(x) } } )
        geneFilterOnStrand = unlist(geneFilterOnStrand)
        return(geneFilterOnStrand)
      })
      myFinalFiltered = unlist(geneFilterOnChr, recursive = TRUE)
      genes = genes[names(genes) != geneId]
      for(newGeneName in names(myFinalFiltered)){
        index = which( names(myFinalFiltered)== newGeneName)
        newGeneId = paste0(geneId,".", index)
        genes[[newGeneId]] = myFinalFiltered[[newGeneName]]
      }
    }
  }
  return(genes)
}    




###################################################################################################
#################################### END FONCTIONS  ###############################################
###################################################################################################

if(!file.exists("./Rsave/introns.rds")){
  dir.create("./Rsave/")
  print("Objects Construction...")
#On initialise ensembl a NULL si jamais aucun dataSet n'est donné pour biomaRt
 ensembl = NULL;

#On importe le GTF
 input = import(gtf);
 #On récupère les annotation BiomaRt SI un dataSet bioMart a été indiquer
 
# input = input[1:2000,];

 if(!is.null(opt$biomaRtDataSet)){
      ensembl = useMart("ensembl",dataset=as.character(opt$biomaRtDataSet));
 }

 #On récupère les lignes dont l'annotation est complète (avec le gene_id)
 annotationGTF = input[!is.na(input$gene_id), ]
 #On récupère que les exons:
 annotationGTF = annotationGTF[annotationGTF$type == "exon", ]
 print("Exons Construction ....")
 #########################################  EXONS ##################################################
 exons = annotationGTF;
 
 ###############################  EXONS BY GENES  ########################################################
 print("ExonsByGenes Construction ...")
 exonsByGenes = split(annotationGTF, annotationGTF$gene_id);
 exonsByGenes = separateTxFromStAndChr(exonsByGenes)
 
 exonsByGenes = unlist(exonsByGenes)
 exonsByGenes$gene_id = names(exonsByGenes)
 exonsByGenes = split(exonsByGenes, exonsByGenes$gene_id)
 
 exonsByGenes = GRangesList(exonsByGenes);
 
 #############################  GENES ##############################################################
 print("Gene Construction ...")
 
 genes = exonsByGenes

 #On utilise lapply pour appliquer sur tout les éléments de la liste afin de créer l'objet GRanges associés
 genes = lapply(X = genes, FUN = function(x){return(GRanges(seqnames = runValue(seqnames(x)), ranges = IRanges(start= min(start(ranges(x))),end = max(end(ranges(x)))), strand = runValue(strand(x))))})
 genes = GRangesList(genes)
 genes = unlist(genes)
 genes$gene_id = names(genes)
 
 print("BiomaRt Construction ...")
 
 if(!is.null(opt$biomaRtDataSet)){
 
 #GENE SYMBOL
     myDataGeneID = as.data.frame(sapply(strsplit(x = genes$gene_id, split = "[.]"), '[', 1));
     colnames(myDataGeneID) = "entrezgene";
     #biomaRt
     biomartGeneSymbol = biomaRt::getBM(attributes=c('entrezgene', 'hgnc_symbol'), filters ='entrezgene', values = myDataGeneID$entrezgene, mart = ensembl) 
     biomartGeneSymbol <- data.frame(lapply(biomartGeneSymbol, as.character), stringsAsFactors = TRUE);
     colnames(biomartGeneSymbol) = c("entrezgene", "hgnc_symbol")
     #filterResults car biomaRt peut renvoyer plusieurs résultata pour un même id, on prend le premier
     biomartGeneSymbol = biomartGeneSymbol[!duplicated(biomartGeneSymbol$entrezgene), ]
     output = dplyr::left_join(myDataGeneID, as.data.frame(biomartGeneSymbol), by = "entrezgene")
     #On associe les identifiants
     genes$gene_symbol = output$hgnc_symbol;
 
 #GENE BIOTYPE
     biomartGeneBiotype = biomaRt::getBM(attributes=c('entrezgene', 'gene_biotype'), filters ='entrezgene', values = myDataGeneID$entrezgene, mart = ensembl) 
     biomartGeneBiotype <- data.frame(lapply(biomartGeneBiotype, as.character), stringsAsFactors = TRUE);
     colnames(biomartGeneBiotype) = c("entrezgene", "gene_biotype");
     #filterResults car biomaRt peut renvoyer plusieurs résultata pour un même id, on prend le premier
     biomartGeneBiotype = biomartGeneBiotype[!duplicated(biomartGeneBiotype$entrezgene), ]
     output <-  dplyr::left_join(myDataGeneID, as.data.frame(biomartGeneBiotype), by = "entrezgene")
     #On associe les identifiants
     genes$gene_biotype = output$gene_biotype;
 }else{
   genes$gene_symbol = "";
   genes$gene_biotype = ""; 
 }
 
 ###############################  TRANSCRITS  ###########################################################
 
 print("Transcrit Construction ...")
 sourceTranscrits = unlist(exonsByGenes, use.names = FALSE)
 
 transcripts = split(sourceTranscrits, sourceTranscrits$transcript_id)

 #On récupère les gene_id ici car je n'arrive pas a les mettre direct dans le GRange
 geneIDForTranscripts = lapply(transcripts, FUN = function(x){return(unique(x$gene_id))})
 transcripts = lapply(X = transcripts, FUN = function(x){return(GRanges(seqnames = runValue(seqnames(x)), ranges = IRanges(start= min(start(ranges(x))),end = max(end(ranges(x)))), strand = runValue(strand(x))))})
 transcripts = GRangesList(transcripts)
 transcripts = unlist(transcripts)
 transcripts$transcript_id = names(transcripts)
 transcripts$gene_id = geneIDForTranscripts
 
 
 

 
 ###############################  EXONS BY TRANSCRITS  ###################################################
 print("exonsByTranscripts Construction ...")
 
 exonsByTranscripts = split(sourceTranscrits, sourceTranscrits$transcript_id)
 exonsByTranscripts = GRangesList(exonsByTranscripts);

 genesWithoutStrand = genes[,NULL]
 strand(genesWithoutStrand) = "*"

 introns = getIntrons(exonsByGenes, genes);
 
 
 saveRDS(object = exons, file = "./Rsave/exons.rds")
 saveRDS(object = genes, file = "./Rsave/genes.rds")
 saveRDS(object = transcripts, file = "./Rsave/transcripts.rds")
 saveRDS(object = exonsByGenes, file = "./Rsave/exonsByGenes.rds")
 saveRDS(object = exonsByTranscripts, file = "./Rsave/exonsByTranscripts.rds")
 saveRDS(object = genesWithoutStrand, file = "./Rsave/genesWithoutStrand.rds")
 saveRDS(object = introns, file = "./Rsave/introns.rds")
 

}else{
 
################################### GET_OBJECTS  ##############################################################

 
 exons = readRDS("./Rsave/exons.rds");
 genes = readRDS("./Rsave/genes.rds");
 transcripts = readRDS("./Rsave/transcripts.rds");
 exonsByGenes = readRDS("./Rsave/exonsByGenes.rds");
 exonsByTranscripts = readRDS("./Rsave/exonsByTranscripts.rds");
 genesWithoutStrand = readRDS("./Rsave/genesWithoutStrand.rds");
 introns = readRDS("./Rsave/introns.rds");


 #####################################  ANNOTATION  ######################################################

} 
 
print("Begin Annotation of candidates ....") 
confidenceWindow = 0; 

data = read.table(mergeTable, header = TRUE, sep = "\t", comment.char = "");

AllTranscriptWithoutStrandWithIds = transcripts
strand(AllTranscriptWithoutStrandWithIds) = "*"

class = TranscriptIds = GeneId = GeneSymbol = GeneBioType = ElementStart = ElementEnd = StartRank = EndRank = ExactMatch = ExactMatchIds = c();
 # On parcours chaque ligne de fichier dans l'ordre et on ajoute les annotations dans l'ordre
 for(i in 1:nrow(data)){
  # print(i)
  myAnnotation = annotate(data[i,]$Start, data[i,]$End, as.character(data[i,]$Chromosome), as.character(data[i,]$Strand));
  class = c(class, myAnnotation[["Class"]])
  TranscriptIds = c(TranscriptIds, myAnnotation[["TranscriptIds"]]);
  GeneId = c(GeneId, concatRes(myAnnotation[["GeneId"]]));
  GeneSymbol = c(GeneSymbol, concatRes(myAnnotation[["GeneSymbol"]]));
  GeneBioType = c(GeneBioType, concatRes(myAnnotation[["GeneBioType"]]));
  ElementStart = c(ElementStart, concatRes(myAnnotation[["ElementStart"]]));
  ElementEnd = c(ElementEnd, concatRes(myAnnotation[["ElementEnd"]]));
  StartRank = c(StartRank, concatRes(myAnnotation[["RankStart"]]));
  EndRank = c(EndRank, concatRes(myAnnotation[["RankEnd"]]));
  ExactMatch = c(ExactMatch, myAnnotation[["ExactMatch"]]);
  ExactMatchIds = c(ExactMatchIds, concatRes(myAnnotation[["ExactMatchTxIds"]]));
}

data$Class = class;
data$TranscriptIds = TranscriptIds;
data$GeneId = GeneId;
data$GeneSymbol = GeneSymbol;
data$GeneBioType = GeneBioType;
data$ElementStart = ElementStart;
data$ElementEnd = ElementEnd;
data$StartRank = StartRank;
data$EndRank = EndRank;
data$ExactMatch = ExactMatch;
data$ExactMatchIds = ExactMatchIds;

write.table(data, file = paste0(path, condition, "_Annotate.csv"),row.names=FALSE, na="", sep = "\t", dec = ".", quote = FALSE)





