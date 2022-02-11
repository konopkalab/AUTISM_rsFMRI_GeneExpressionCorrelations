# Covariate plot
rm(list=ls())
suppressPackageStartupMessages({
library(here)
library(tidyverse)
library(data.table)
})

dir.create("supp_tables")

tab1 <- read.table("integrative_results/ASD_fMRI_Genes.txt",header=T,sep="\t")
tab2 <- read.table("integrative_results/ASD_fMRI_DevClustered.txt",header=T,sep="\t")

tab3 <- tab1 %>%
        filter(Direction == "TopRsq") %>%
        rename(TopRsq = Direction)

tab4 <- tab1 %>%
        filter(!Direction %in% c("TopRsq","AllGene"))

tmp <- list(tab2,tab3,tab4)
df <- Reduce(function(x, y) {
    merge(x, y, all=TRUE, by="Gene")
}, tmp)

df <- full_join(tab2,tab3,by="Gene") %>%
        filter(!Direction == "AllGene") %>%
        arrange(Direction)

wgcna <- list(WGCNA = df)


Neuron <- read.table("rawdata/fALFF_genes_Pearson.txt",header=T)
Neuron <- data.frame(Gene=Neuron$fALFF_genes, WangEtAl=rep("WangEtAl",nrow(Neuron)))
Rich <- read.table("rawdata/Richiardi_science.txt",header=T)
Rich <- data.frame(Gene=Rich$Gene, RichiardiEtAl=rep("RichiardiEtAl",nrow(Rich)))

fMRI=list(WangEtAl = Neuron, RichiardiEtAl = Rich)

# Load data for Data Base
psydge <- load(here("rawdata","geneset_for_enrichment","PsychENCODE_DEGs.RData")) %>%
            get()
new_names <- names(psydge)
psydge <- map2(psydge, new_names, ~setnames(.x, 'Class', .y))

psymod <- load(here("rawdata","geneset_for_enrichment", "PsychEncode_Modules.RData")) %>%
            get()
new_names <- names(psymod)
psymod <- map2(psymod, new_names, ~setnames(.x, 'Mod', .y))

sfari <- load(here("rawdata","geneset_for_enrichment","ASD_SFARI.RData")) %>%
            get()
new_names <- names(sfari)
sfari <- map2(sfari, new_names, ~setnames(.x, 'Class', .y))

scBA <- load(here("rawdata","geneset_for_enrichment", "Allen_MultiReg_CellMarkers_GeneSet.RData")) %>%
            get()
new_names <- names(scBA)
scBA <- map2(scBA, new_names, ~setnames(.x, 'cluster', .y))

reg <- load(here("rawdata","geneset_for_enrichment", "DEGs_ByReg.RData")) %>%
            get()
new_names <- names(reg)
reg <- map2(reg, new_names, ~setnames(.x, 'ID', .y))

l <- c(wgcna,fMRI,sfari,reg,psydge,psymod,scBA)
res <- Reduce(function(x, y) {
    merge(x, y, all=TRUE, by="Gene")
}, l)

database=res[!(is.na(res$Direction)),]
database <- database[!(duplicated(database)),]
openxlsx::write.xlsx(database, file = "supp_tables/Table_S2.xlsx", 
                    colNames = TRUE, borders = "columns",overwrite=TRUE)



