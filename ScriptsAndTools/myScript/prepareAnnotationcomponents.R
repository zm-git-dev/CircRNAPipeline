library("optparse");
library("magrittr");
library("rtracklayer");
library("biomaRt");
library("GenomicRanges")






getArguements = function(){
  option_list = list(
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
gtf = as.character(opt$gtf);
idBDD = as.character(opt$idBDD);

# biomaRtDataSet est directement utilisé dans le if pour savoir s'il a été envoyé sinon il n'y aura pas de geneSymbol et GeneBioType

###################################################################################################
################################### START FONCTIONS  ##############################################
###################################################################################################





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
  
  ###############################  EXONS BY GENES  ##################################################
  print("ExonsByGenes Construction ...")
  exonsByGenes = split(annotationGTF, annotationGTF$gene_id);
  exonsByGenes = separateTxFromStAndChr(exonsByGenes)
  
  exonsByGenes = unlist(exonsByGenes)
  exonsByGenes$gene_id = names(exonsByGenes)
  exonsByGenes = split(exonsByGenes, exonsByGenes$gene_id)
  
  exonsByGenes = GRangesList(exonsByGenes);
  
  # Separation des Locus par transcripts
  ####### On détermine les transcripts une première fois  ####### 
  PresourceTranscrits = unlist(exonsByGenes, use.names = FALSE)
  Pretranscripts = split(PresourceTranscrits, PresourceTranscrits$transcript_id)
  
  #On récupère les gene_id ici car je n'arrive pas a les mettre direct dans le GRange
  PregeneIDForTranscripts = lapply(Pretranscripts, FUN = function(x){return(unique(x$gene_id))})
  Pretranscripts = lapply(X = Pretranscripts, FUN = function(x){return(GRanges(seqnames = runValue(seqnames(x)), ranges = IRanges(start= min(start(ranges(x))),end = max(end(ranges(x)))), strand = runValue(strand(x))))})
  Pretranscripts = GRangesList(Pretranscripts)
  Pretranscripts = unlist(Pretranscripts)
  Pretranscripts$transcript_id = names(Pretranscripts)
  Pretranscripts$gene_id = as.character(PregeneIDForTranscripts)
  
  # On découpe alors par gènes
  transcriptsByGene = split(Pretranscripts, Pretranscripts$gene_id)
  # On cherche les transcirpt d'un même gene_id qui ne se chevauche pas
  cc = lapply(names(transcriptsByGene), FUN = function(x){
    length(transcriptsByGene[[x]] %>% reduce())
    })
  # On détermines les gènes concernés
  GenesOfInterest = names(transcriptsByGene)[unlist(cc) > 1]
  for(GenesOfInterest_id in GenesOfInterest){
    #On récupère tout les exons de ce gènes
    AllGeneExons = exonsByGenes[[GenesOfInterest_id]]
    # On détermine les deux future régions du gène c a d la on les transcripts ne s'overlap pas
    reducedTX = transcriptsByGene[[GenesOfInterest_id]] %>% reduce()
    # pour chacune de ces régions
    for(i in 1:(length(reducedTX))){
      # On récupère les exons concernés dans cette région
      exonsSelect = AllGeneExons[subjectHits(findOverlaps(reducedTX[i], AllGeneExons)), ]
      # On leur attribut un nouveau gene_id
      newGeneId = paste0(c(exonsSelect[1,]$gene_id, i), collapse = ".")
      exonsSelect$gene_id = newGeneId
      # On l'ins_re dans la liste
      exonsByGenes[[newGeneId]] = exonsSelect
    }
  }
  # On suprime les anciens
  exonsByGenes = exonsByGenes[! names(exonsByGenes) %in%  GenesOfInterest]
  
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
  
}
