suppressPackageStartupMessages({
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(UpSetR)
library(RColorBrewer)
library(tidyverse)
library(VennDiagram)
library(ggVennDiagram)
library(clusterProfiler)
library(RRHO2)
library(here)
})



# gene Onto
genes <- read.table(here("integrative_results","ASD_fMRI_Genes.txt"),header=T,sep="\t")
universe <- read.table("UTILS/BrainExpressed_Background.txt",header=T)

all <- genes %>%
        filter(Direction == "TopRsq")

GOI_Uni <- bitr(as.character(universe$x),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db)


GOI <- bitr(as.character(all$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db)

GeneOnto_TopRsq <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID",
                     universe =  unique(GOI_Uni$ENTREZID),
                     OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
                     ont = "MF", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 1, 
                     qvalueCutoff = 1,  
                     readable = TRUE) %>%
					 filter(Count > 5)

pdf("revision/GeneOnto_TopRsq.pdf",width=6,height=5,useDingbats=FALSE)
clusterProfiler::dotplot(GeneOnto_TopRsq, showCategory=4, color="pvalue")
dev.off()


openxlsx::write.xlsx(as.data.frame(GeneOnto_TopRsq), 
                     file = "revision/GeneOnto_TopRsq.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto_TopRsq", 
                     overwrite = TRUE)

all <- genes %>%
        filter(Direction == "AllGene")


GOI <- bitr(as.character(all$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db)

GeneOnto_Adult <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     #universe =  unique(GOI_Uni$ENTREZID),
                     OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
                     ont = "MF", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 1, 
                     qvalueCutoff = 1, 
                     readable = TRUE)%>%
					 filter(Count > 5)

pdf("revision/GeneOnto_AllGenes.pdf",width=10,height=5,useDingbats=FALSE)
dotplot(GeneOnto_Adult, showCategory=4, color="pvalue")
dev.off()


openxlsx::write.xlsx(as.data.frame(GeneOnto_Adult), 
                     file = "revision/GeneOnto_AllGenes.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto_AllGenes", 
                     overwrite = TRUE)
