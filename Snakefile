configfile: "config.json"
include: "rules/compute_mergedReads.rules"
include: "rules/bowtie2_index.rules" 
include: "rules/bowtie2_align.rules"
include: "rules/compteUnmappedAnchors.rules"
include: "rules/computeFindCirc.rules"
include: "rules/star_index.rules"
include: "rules/star_aln_with_chimeric.rules"
include: "rules/parse_chimeric_junction.rules"
include: "rules/eliminateArtefact.rules"
include: "rules/mergePrediction.rules"
include: "rules/annotate_circRna.rules"
include: "rules/compute_circRNA_finder.rules"
include: "rules/createRDSFiles.rules"
include: "rules/catchProcessInformations.rules"
include: "rules/createRPackageForGenome.rules"
include: "rules/runHtSeqCount.rules"
include: "rules/create_annotation_components.rules"
include: "rules/launchtest.rules"
include: "rules/prepare_unmerge_reads.rules"
include: "rules/removeGenomeInSharedMemory.rules"
localrules: compute_mergedReads

rule all:
	input:
		"dataAnalyseCirc.rds",
		"dataAnalyseProcess.rds",
		"dataHTSEQcount.rds"
