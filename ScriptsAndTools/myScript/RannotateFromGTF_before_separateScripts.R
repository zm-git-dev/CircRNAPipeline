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
                help="Name of the biomaRt dataset to use for annotation", metavar="character"),
    make_option(c("--idBDD"), type = "character", default = NULL,
                help="Name source of gene_id RefSeq or Ensembl", metavar="character")
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
idBDD = as.character(opt$idBDD);

# biomaRtDataSet est directement utilisé dans le if pour savoir s'il a été envoyé sinon il n'y aura pas de geneSymbol et GeneBioType

#### Avec RtrackLayer

# gtf = "./newGTF.gtf";
# path = "./Janvier2014/";
# condition ="Senescence";
# biomaRtDataset = "hsapiens_gene_ensembl";
# mergeTable = "./Janvier2014/Senescence/SenescenceFinalMergeTest.csv";
 
# idBDD = "Ensembl"

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
  V_ElementLocus = V_ELementTx = V_ExactMatch = V_ExactMatchTxIds = V_Class = c()
      # Pour chaque gène : DANS L'ORDRE
      for(i in (1:length(Genes))){
          gene = Genes[i,]
          # On cherche les élements génomique en Strand et Un_Strand
          PreviousV_Element_Strand = computeAnnotation(GR, gene$gene_id, exonsByGenes, intronsInGenes, "@E", "@I", as.character(strand(gene)), strand, 0)
          PreviousV_Element_Unstrand =   computeAnnotation(GR, gene$gene_id, exonsByGenes, intronsInGenes, "@AE", "@AI", as.character(strand(gene)), ifelse(strand == "+", "-", "+"), 10)
          # En fait pour le gène i au strat ou au end  si il y a un élément génomique il n' en a qu'un donc on regarde si la concat est > 0 si non c'est donc que rien n'as mappé on met NA
          PreviousV_Element = c(PreviousV_Element_Strand$Element, PreviousV_Element_Unstrand$Element)
          V_ElementLocus = c(V_ElementLocus, ifelse(length(PreviousV_Element) > 0, yes = PreviousV_Element, no = "NA"))
          
          V_Class = c(V_Class, (PreviousV_Element_Strand$Class + PreviousV_Element_Unstrand$Class))
          
          TxForGene_i = Tx[Tx$gene_id == gene$gene_id,]          
          if(length(TxForGene_i) > 0){
              for(j in (1:length(TxForGene_i))){
                transcript = TxForGene_i[j, ]
                # On fait pareil pour tout les transcripts associés au gène
                PreviousV_Element = c(computeAnnotation(GR, transcript$transcript_id, exonsByTranscripts, intronsInTranscripts, "@E", "@I", as.character(strand(transcript)), strand, 0)$Element,
                                      computeAnnotation(GR, transcript$transcript_id, exonsByTranscripts, intronsInTranscripts, "@AE", "@AI", as.character(strand(transcript)), ifelse(strand == "+", "-", "+"), 10)$Element)
                
                V_ELementTx = c(V_ELementTx, ifelse(length(PreviousV_Element) > 0, yes = PreviousV_Element, no = "NA"))
                
              }
          }
          # On fait le test de l'exact match à l'échelle du gène surlequel on travaille
          TxKeepedExactMatchBool = testIfExactMatch(TxForGene_i, GR, test)
          V_ExactMatch = c(V_ExactMatch, ifelse(sum(TxKeepedExactMatchBool) != 0, TRUE, FALSE))
          V_ExactMatchTxIds = c(V_ExactMatchTxIds, TxForGene_i[TxKeepedExactMatchBool, ]$transcript_id)
      }
      V_Class = as.vector(unlist(TableOfClass[as.character(V_Class)]))
      
  resultat[[paste0("ElementLocus",test)]] = V_ElementLocus
  resultat[[paste0("ElementTranscripts",test)]] = V_ELementTx
  resultat[[paste0("Class",test)]] = V_Class
  resultat[[paste0("ExactMatch",test)]] = V_ExactMatch
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
  
  AllGenesM = base::union(genes[subjectHits(findOverlaps(GRStartWithoutStrand, genes)),]$gene_id, genes[subjectHits(findOverlaps(GREndWithoutStrand, genes)),]$gene_id)
  resultat = newAnnotateGR(GRStartWithoutStrand, "Start", strand, resultat, AllGenesM);
  resultat = newAnnotateGR(GREndWithoutStrand, "End", strand, resultat, AllGenesM);
  resultat[is.na(resultat)] = ""
  Genes = genes[genes$gene_id %in% AllGenesM, ]
  
  nbGenes = ifelse(length(AllGenesM) > 0, length(AllGenesM), 1)
  
  newDataAnnotate = dataLine[rep(1, nbGenes),]
  newDataAnnotate$Class = resumeAnnotationFactor(resultat, "Class", "-")
  newDataAnnotate$TranscriptIds = resumeAllAnnotationFactor(resultat, "TranscriptIds", ",")
  newDataAnnotate$GeneId = paste0(Genes$gene_id, "")
  newDataAnnotate$GeneSymbol = paste0(Genes$gene_symbol, "")
  newDataAnnotate$GeneBioType = paste0(Genes$gene_biotype, "")
  newDataAnnotate$ElementLocusStart = resultat[["ElementLocusStart"]]
  newDataAnnotate$ElementLocusEnd = resultat[["ElementLocusEnd"]]
  newDataAnnotate$ElementTranscriptsStart = concatRes(resultat[["ElementTranscriptsStart"]], ",")
  newDataAnnotate$ElementTranscriptsEnd = concatRes(resultat[["ElementTranscriptsEnd"]], ",")
  newDataAnnotate$ExactMatch =  resultat[["ExactMatchStart"]] & resultat[["ExactMatchEnd"]]
  newDataAnnotate$ExactMatchIds = resumeAllAnnotationFactor(resultat, "ExactMatchTxIds", ",")

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
  resultat[[paste0("TranscriptIds", test)]] = resultat[[paste0("GeneId",test)]] = resultat[[paste0("GeneSymbol",test)]] = resultat[[paste0("GeneBioType",test)]] = resultat[[paste0("ElementLocus",test)]] = resultat[[paste0("ElementTranscripts",test)]] = resultat[[paste0("ElementLocus", test)]] = resultat[[paste0("ElementTranscripts", test)]] = resultat[[paste0("ExactMatchTxIds", test)]] = "";
  resultat[[paste0("ExactMatch", test)]] = FALSE
  return(resultat);
}


# Crée un objet GRanges représentant les introns à paritr des coordonnées des exons de chaque Transcripts
getIntronsInTranscripts  = function(exonsByTranscripts, transcripts){
  GRIntrons = GRanges();
  for(txName in names(exonsByTranscripts)){
    exons <- exonsByTranscripts[[txName]] %>% sort();
    mcols(exons) = NULL
    if(length(exons) > 1){
      currentGRTx = GRanges(seqnames=seqnames(exons)[1], ranges=IRanges(start=min(start(exons)), end=max(end(exons))), strand=strand(exons)[1])  
      db = disjoin(c(exons, currentGRTx));
      ints = db[countOverlaps(db, exons) == 0]
      #Add an ID
      if(length(ints) > 0){
        if(as.character(strand(currentGRTx)) == "-") {
          ints$intron_id = c(length(ints):1)
        } else {
          ints$intron_id = c(1:length(ints))
        }
        
        ints$transcript_id = txName
        ints$gene_id = as.character(transcripts[transcripts$transcript_id == txName,]$gene_id)
        GRIntrons = append(GRIntrons, ints);
      }  
    }
  }
  return(GRIntrons);
}



# Crée un objet GRanges représentant les introns à paritr des coordonnées des exons de chaque Locus
getIntronsInGenes  = function(exonsByGenes, genes){
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

#On importe le GTF
input = import(gtf);

if(!file.exists("./Rsave/intronsInGenes.rds")){
  dir.create("./Rsave/")
  print("Objects Construction...")
#On initialise ensembl a NULL si jamais aucun dataSet n'est donné pour biomaRt
 ensembl = NULL;

# input = import(gtf);


#On récupère les annotation BiomaRt SI un dataSet bioMart a été indiquer
 
 # On détermine a quel type va traiter si il y a requête biomaRt
 IdType = list(RefSeq = "entrezgene", Ensembl = "ensembl_gene_id")
 
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
 
 if(!is.null(opt$biomaRtDataSet)){
   
   print("BiomaRt Construction ...")
   # On récupère le nom du champ associé au type d'ID dans le GTF
   idToRequest = IdType[[idBDD]]
   
 #GENE SYMBOL
     myDataGeneID = as.data.frame(sapply(strsplit(x = genes$gene_id, split = "[.]"), '[', 1));
     colnames(myDataGeneID) = "geneID";
     #biomaRt
     biomartGeneSymbol = biomaRt::getBM(attributes=c(idToRequest, 'hgnc_symbol'), filters = idToRequest, values = myDataGeneID$geneID, mart = ensembl) 
     biomartGeneSymbol <- data.frame(lapply(biomartGeneSymbol, as.character), stringsAsFactors = TRUE);
     colnames(biomartGeneSymbol) = c("geneID", "hgnc_symbol")
     #filterResults car biomaRt peut renvoyer plusieurs résultata pour un même id, on prend le premier
     biomartGeneSymbol = biomartGeneSymbol[!duplicated(biomartGeneSymbol$geneID), ]
     output = dplyr::left_join(myDataGeneID, as.data.frame(biomartGeneSymbol), by = "geneID")
     #On associe les identifiants
     genes$gene_symbol = output$hgnc_symbol;
 
 #GENE BIOTYPE
     biomartGeneBiotype = biomaRt::getBM(attributes=c(idToRequest, 'gene_biotype'), filters =idToRequest, values = myDataGeneID$geneID, mart = ensembl) 
     biomartGeneBiotype <- data.frame(lapply(biomartGeneBiotype, as.character), stringsAsFactors = TRUE);
     colnames(biomartGeneBiotype) = c("geneID", "gene_biotype");
     #filterResults car biomaRt peut renvoyer plusieurs résultata pour un même id, on prend le premier
     biomartGeneBiotype = biomartGeneBiotype[!duplicated(biomartGeneBiotype$geneID), ]
     output <-  dplyr::left_join(myDataGeneID, as.data.frame(biomartGeneBiotype), by = "geneID")
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

 # Pour les introns, on prend le même format que pour les exons
 intronsInTranscripts = getIntronsInTranscripts(exonsByTranscripts, transcripts);
 intronsInGenes = getIntronsInGenes(exonsByGenes, genes);

 
 saveRDS(object = exons, file = "./Rsave/exons.rds")
 saveRDS(object = genes, file = "./Rsave/genes.rds")
 saveRDS(object = transcripts, file = "./Rsave/transcripts.rds")
 saveRDS(object = exonsByGenes, file = "./Rsave/exonsByGenes.rds")
 saveRDS(object = exonsByTranscripts, file = "./Rsave/exonsByTranscripts.rds")
 saveRDS(object = genesWithoutStrand, file = "./Rsave/genesWithoutStrand.rds")
 saveRDS(object = intronsInTranscripts, file = "./Rsave/intronsInTranscripts.rds")
 saveRDS(object = intronsInGenes, file = "./Rsave/intronsInGenes.rds")

}else{
 
################################### GET_OBJECTS  ##############################################################

 
 exons = readRDS("./Rsave/exons.rds");
 genes = readRDS("./Rsave/genes.rds");
 transcripts = readRDS("./Rsave/transcripts.rds");
 exonsByGenes = readRDS("./Rsave/exonsByGenes.rds");
 exonsByTranscripts = readRDS("./Rsave/exonsByTranscripts.rds");
 genesWithoutStrand = readRDS("./Rsave/genesWithoutStrand.rds");
 intronsInTranscripts = readRDS("./Rsave/intronsInTranscripts.rds");
 intronsInGenes = readRDS("./Rsave/intronsInGenes.rds");

 #####################################  ANNOTATION  ######################################################

} 

intronsInTranscripts = split(intronsInTranscripts, intronsInTranscripts$transcript_id); 
intronsInGenes = split(intronsInGenes, intronsInGenes$gene_id);

print("Begin Annotation of candidates ....") 
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


