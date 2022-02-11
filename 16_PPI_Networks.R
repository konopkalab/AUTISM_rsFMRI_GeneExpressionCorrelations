rm(list=ls())
suppressPackageStartupMessages({
library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(here)
library(STRINGdb)
library(igraph)
library(rio)
library(effsize)
library(httr)
})
source("UTILS/Utils.R")

#set_config(
#  use_proxy(url="proxy.swmed.edu", port=3128,username="S157784",password="Arianna2103!")
#)

# Create Output directories
folder_names <- c("networking")
sapply(folder_names, dir.create)

#Load the final gene list of differentially correlated genes
finalgenes = read.table("integrative_results/ASD_fMRI_DevClustered.txt", header = T)
names <- as.character(unique(finalgenes$Direction))

# Create string db
string_db <- STRINGdb$new( version="11", species=9606, score_threshold=400, input_directory="" )

bkg <- read.table("UTILS/BrainExpressed_Background.txt")
string_db$set_background(as.character(bkg$x))

# Loop across all the cool modules
input <- list()
mapped <- list()
hits <- list()
enrich <- list()
pdfname <- list()
for(i in 1:length(names))
	{
		pdfname[[i]]<-paste0("networking/",names[[i]],"_PPI_network_BrainExp.pdf") 
		pdf(file=pdfname[[i]],,width = 5,height=5)
		input[[i]] <- finalgenes %>% filter(Direction==names[[i]])
		mapped[[i]] <- string_db$map(input[[i]], "Gene", removeUnmappedRows = TRUE ) #%>% as.data.frame()
		enrich[[i]] <- string_db$get_ppi_enrichment(mapped[[i]]$STRING_id[1:nrow(input[[i]])])
 		hits[[i]] <- mapped[[i]]$STRING_id[1:nrow(input[[i]])]
		string_db$plot_network(hits[[i]])
		dev.off()
	}


# Coexnet ppi
finalgenes = read.table("integrative_results/ASD_fMRI_DevClustered.txt", header = T)

adult <- finalgenes %>% filter(Direction == "Adult")

ppi <- ppiNet(molecularIDs = as.character(adult$Gene),
                        Score = 400,
				evidence = c("neighborhood","coexpression","experiments","combined_score"))


pdf("networking/Adult_IGRAPH.pdf",width=4,height=4)
ppi <- delete.vertices(ppi,which(degree(ppi)<2))
layoutFR <- layout_with_fr(ppi,maxiter = 500)
plot.igraph(ppi,
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=3,
            vertex.label.color="black",
            vertex.color="royalblue",
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=adjustcolor("grey", alpha.f = .5))
dev.off()


adult <- finalgenes %>% filter(Direction == "EarlyDev")

ppi <- ppiNet(molecularIDs = as.character(adult$Gene),
				Score = 400,
				evidence = c("neighborhood","coexpression","experiments","combined_score"))

pdf("networking/EarlyDev_IGRAPH.pdf",width=4,height=4)
ppi <- delete.vertices(ppi,which(degree(ppi)<2))
layoutFR <- layout_with_fr(ppi,maxiter = 500)
plot.igraph(ppi,
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=3,
            vertex.label.color="black",
            vertex.color="royalblue",
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=adjustcolor("grey", alpha.f = .5))
dev.off()


adult <- finalgenes %>% filter(Direction == "Stable")

ppi <- ppiNet(molecularIDs = as.character(adult$Gene),
				Score = 400,
				evidence = c("neighborhood","coexpression","experiments","combined_score"))

pdf("networking/Stable_IGRAPH.pdf",width=4,height=4)
ppi <- delete.vertices(ppi,which(degree(ppi)<2))
layoutFR <- layout_with_fr(ppi,maxiter = 500)
plot.igraph(ppi,
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=3,
            vertex.label.color="black",
            vertex.color="royalblue",
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=adjustcolor("grey", alpha.f = .5))
dev.off()

