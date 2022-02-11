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

load("integrative_results/DiffCor_Significant.RData")

# BarPlot

pdf(file="revision/Barplot_DCgene_SignCor.pdf",width=4,height=4)
DiffCor_Sign_Final %>% 
group_by(Direction) %>% 
tally() %>%
 ggbarplot("Direction", "n",
   fill = "Direction", color = "Direction",
   palette = c("#FDE725FF", "#440154FF"),
   label = TRUE, lab.vjust = -1.6) +
theme(legend.position="none") +
          rotate_x_text(angle = 45)+
          xlab("")+
          ylab("")
dev.off()


# Scatter Plots
pdf(file="revision/ScatterPlot_CTL_Rho.pdf",width=4,height=4)
DiffCor_Sign_Final %>% 
ggscatter( 
            x = "Rho_CTL_ReHo", 
            y = "Rho_CTL_fALFF",
            color = "Direction",
            palette=c("#FDE725FF", "#440154FF"),
            size = 1,
            alpha=0.3,
            shape=19,
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = -0.4,lable.y=0.2, label.sep = "\n"))+
      xlab("ReHo's rho")+ 
      ylab("fALFF's rho")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      theme(legend.position="none")+
      ylim(-0.4,0.4)+
      xlim(-0.4,0.4)
dev.off()

pdf(file="revision/ScatterPlot_ASD_Rho.pdf",width=4,height=4)
DiffCor_Sign_Final %>% 
ggscatter(
            x = "Rho_ASD_ReHo", 
            y = "Rho_ASD_fALFF",
            color = "Direction",
            palette=c("#FDE725FF", "#440154FF"),
            size = 1,
            alpha=0.3,
            shape=19,
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = -0.4,lable.y=0.2, label.sep = "\n"))+
      xlab("ReHo's rho")+ 
      ylab("fALFF's rho")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      theme(legend.position="none")+
      ylim(-0.4,0.4)+
      xlim(-0.4,0.4)
dev.off()

pdf(file="revision/ScatterPlot_DiffCorZ.pdf",width=4,height=4)
p <- CTL_Sign_Final %>%
         mutate(Sign = if_else(Gene %in% DiffCor_Sign_Final$Gene, "Sign", "NotSign"), 
                MeanZ = (DiffCor_ReHo_Z+DiffCor_fALFF_Z)/2) %>% 
         ggscatter(
            x = "DiffCor_fALFF_Z", 
            y = "DiffCor_ReHo_Z",
            color = "Sign",
            palette=c("grey","black"),
            size = "MeanZ",
            alpha=0.3,
            shape=19,
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 0.1,label.y=4, label.sep = "\n"))+
      xlab("ReHo's Diff Cor Z")+ 
      ylab("fALFF's Diff Cor Z") +
      theme(legend.position="none") +
      scale_size(range = c(0, 3))
ggMarginal(p,type = "density")
dev.off()

# Vulcano Plot
pdf(file="revision/Volcano_DiffCorZ.pdf",width=8,height=8)

top_labelled <- CTL_Sign_Final %>% 
                  mutate(Sign = if_else(Gene %in% DiffCor_Sign_Final$Gene, "Sign", "NotSign"), 
                MeanZ = (DiffCor_ReHo_Z+DiffCor_fALFF_Z)/2,
                LOG = -log10(DiffCor_Comb_P)) %>%
                  na.omit() %>%
                  arrange(MeanZ) %>%
                  top_n(n = 15, wt = MeanZ)

CTL_Sign_Final %>%
         mutate(Sign = if_else(Gene %in% DiffCor_Sign_Final$Gene, "Sign", "NotSign"), 
                MeanZ = (DiffCor_ReHo_Z+DiffCor_fALFF_Z)/2,
                LOG = -log10(DiffCor_Comb_P)) %>%
               ggscatter(
                           x = "MeanZ", 
                           y = "LOG",
                           color = "Sign",
                           palette=c("grey","red"),
                           size = 1,
                           alpha=0.3,
                           shape=19)+
                     xlab("Diff Cor Z")+ 
                     ylab("-log10(FDR)")+
                     geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
                     #geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
                     #geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
                     geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
                     geom_text_repel(data = top_labelled, 
                                     mapping = aes(label = Gene), 
                                     size = 5,
                                     box.padding = unit(0.4, "lines"),
                                     point.padding = unit(0.4, "lines"))+
                     theme(legend.position="none")+
                     xlim(0,5)+
                     ylim(0,5)
dev.off()




# gene Onto
genes <- read.table(here("integrative_results","ASD_fMRI_DevClustered.txt"),header=T,sep="\t")
universe <- read.table("UTILS/BrainExpressed_Background.txt",header=T)

all <- genes %>%
        filter(Direction == "EarlyDev")

GOI_Uni <- bitr(as.character(universe$x),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db)


GOI <- bitr(as.character(all$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db)

GeneOnto_EarlyDev <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID",
                     universe =  unique(GOI_Uni$ENTREZID),
                     OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
                     ont = "MF", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 1, 
                     qvalueCutoff = 1,  
                     readable = TRUE) %>%
					 filter(Count > 5)

pdf("revision/GeneOnto_EarlyDev.pdf",width=10,height=5,useDingbats=FALSE)
clusterProfiler::dotplot(GeneOnto_EarlyDev, showCategory=4, color="pvalue")
dev.off()


openxlsx::write.xlsx(as.data.frame(GeneOnto_EarlyDev), 
                     file = "revision/GeneOnto_EarlyDev.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto_EarlyDev", 
                     overwrite = TRUE)

all <- genes %>%
        filter(Direction == "Adult")


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

pdf("revision/GeneOnto_Adult.pdf",width=10,height=5,useDingbats=FALSE)
dotplot(GeneOnto_Adult, showCategory=4, color="pvalue")
dev.off()


openxlsx::write.xlsx(as.data.frame(GeneOnto_Adult), 
                     file = "revision/GeneOnto_Adult.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto_Adult", 
                     overwrite = TRUE)


all <- genes %>%
        filter(Direction == "Stable")


GOI <- bitr(as.character(all$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db)

GeneOnto_Stable <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
                     universe =  unique(GOI_Uni$ENTREZID),
                     ont = "MF", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 1, 
                     qvalueCutoff = 1, 
                     readable = TRUE) %>%
					 filter(Count > 5)

pdf("revision/GeneOnto_Stable.pdf",width=10,height=5,useDingbats=FALSE)
dotplot(GeneOnto_Stable, showCategory=4, color="pvalue")
dev.off()


openxlsx::write.xlsx(as.data.frame(GeneOnto_Stable), 
                     file = "revision/GeneOnto_Stable.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto_Stable", 
                     overwrite = TRUE)