# Covariate plot
rm(list=ls())
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RcisTarget))
suppressPackageStartupMessages(library(igraph))

# Create Output directories
folder_names <- c("motif_enrichment")
sapply(folder_names, dir.create)

#Load the final gene list of differentially correlated genes
finalgenes = read.table("integrative_results/ASD_fMRI_DevClustered.txt", header = T)
names <- as.character(unique(finalgenes$Clusters))
# Load the motid database
# Download the database from "https://resources.aertslab.org/cistarget/"
motifRankings <- importRankings("rawdata/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
data(motifAnnotations_hgnc)

# Loop across all the cool modules
input <- list()
geneLists <- list()
hits <- list()
enrich <- list()
for(i in 1:length(names))
	{
		input[[i]] <- finalgenes %>% filter(Clusters==names[[i]])
		geneLists[[i]] <- list(geneSetName=as.character(input[[i]]$Gene))
		enrich[[i]] <- cisTarget(geneLists[[i]], motifRankings,motifAnnot=motifAnnotations_hgnc,nCores = 10) %>% addLogo()
		write.table(enrich[[i]],paste0("motif_enrichment/",names[[i]],"_Motif_Enrichment.txt"), sep="\t",quote=F)
	}
