rm(list=ls())
suppressPackageStartupMessages({
library(ggplot2)
library(ggjoy)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(GGally)
library(UpSetR)
library(tidyverse)
library(here)
library(ggforce)
library(VennDiagram)
library(DGCA)
library(data.table)
})

# Create Directory
folder_names <- c("integrative_visualizations")
sapply(folder_names, dir.create)

# loading genomic data and formatting expression and demographic data
load(here("rawdata","Expression_Data.RData"))
demo$Region <- gsub("-","_",demo$Region)
demo$ID <- as.factor(paste(demo$Subject,demo$Region,sep="_"))

demoCTL <- demo %>%
            rownames_to_column("OriginalID") %>%
            filter(Diagnosis == "CTL") %>%
            arrange(Region) %>%
            droplevels() %>%
            column_to_rownames("OriginalID")

demoASD <- demo %>%
            rownames_to_column("OriginalID") %>%
            filter(Diagnosis == "ASD") %>%
            arrange(Region) %>%
            droplevels() %>%
            column_to_rownames("OriginalID")

p_ctl <- expAdj[,colnames(expAdj) %in% rownames(demoCTL)]
p_asd <- expAdj[,colnames(expAdj) %in% rownames(demoASD)]
p_ctl <- p_ctl[,match(rownames(demoCTL),colnames(p_ctl))]
colnames(p_ctl) <- demoCTL$ID
p_asd <- p_asd[,match(rownames(demoASD),colnames(p_asd))]
colnames(p_asd) <- demoASD$ID


# Check the gene correlations between values 
load("integrative_results/DiffCor_Tables.RData")

pdf(file="integrative_visualizations/GeneCorrelations_ScatterCor_1.pdf",width=5,height=5)
ggscatter(DiffCor_Combined, x = "Rho_CTL_fALFF", y = "Rho_CTL_ReHo",
   color = "lightgray", size = 0.5, 
   add = "reg.line",  
   add.params = list(color = "red", fill = "blue"), 
   conf.int = TRUE, 
   cor.coef = TRUE, 
   cor.coeff.args = list(method = "pearson", label.sep = "\n"),
   ellipse = TRUE)+ 
xlab("fALFF ~ Gene Expression Rho") + ylab("ReHo ~ Gene Expression Rho") +
geom_hline(yintercept = 0, linetype="dotted", color = "black", size=1)+ 
geom_vline(xintercept = 0, linetype="dotted", color = "black", size=1)
dev.off()


#################### 
### P-value Dist ###
####################
pdf("integrative_visualizations/Pvalue_Distribution_fALFF.pdf",width=5,height=4,useDingbats=FALSE)
DiffCor_Combined %>% 
  select(Pval_CTL_fALFF,Pval_ASD_fALFF) %>%
  melt() %>%
    ggdensity(x = "value",
       color = "variable", 
       fill = "variable",
       palette = c("black", "red"),
       ggtheme = theme_pubr())+
    labs(title="fALFF",x="p-value", y = "")+
    theme_classic()+
    theme(legend.position = c(0.8, 0.8))+
    xlim(0,1)+
    ylim(0,4)
dev.off()

pdf("integrative_visualizations/Pvalue_Distribution_ReHo.pdf",width=5,height=4,useDingbats=FALSE)
DiffCor_Combined %>% 
  select(Pval_CTL_ReHo,Pval_ASD_ReHo) %>%
  melt() %>%
    ggdensity(x = "value",
       color = "variable", 
       fill = "variable",
       palette = c("black", "red"),
       ggtheme = theme_pubr())+
    labs(title="ReHo",x="p-value", y = "")+
    theme_classic()+
    theme(legend.position = c(0.8, 0.8))+
    xlim(0,1)+
    ylim(0,4)
dev.off()


# Load data for visualizations
load(here("integrative_results","Only_CTL_Significant.RData"))
load(here("integrative_results","DiffCor_Significant.RData"))

l <- list(Only_CTL_fALFF,Only_CTL_ReHo)
l2 <- list(DiffCor_Sign_Final)


new_order <- c("fALFF","ReHo")


tmp <- data.frame(ID = c("fALFF","ReHo"),CTL=sapply(l,nrow),DiffCor = sapply(l2,nrow)) %>%
        melt() %>%
        arrange(match(ID, new_order))%>%
        mutate(ID = factor(ID, levels = new_order))

pdf(file="integrative_visualizations/Barplot_NumberOfGenes.pdf",width=4,height=3)
ggbarplot(tmp, 
          "ID", 
          "value",
          fill = "variable", 
          color = "variable", 
          palette = c("black","red"),
          label = TRUE,
          position = position_dodge(0.75))+
          theme(legend.position="right")+
          rotate_x_text(angle = 45)+
          xlab("")+
          ylab("# of Genes")+
          ylim(0,5000)
dev.off()


# File for enrichments
df <- DiffCor_Sign_Final %>% 
            select(Gene, Direction) %>%
            arrange(Direction)

df2 <- DiffCor_Sign_Final %>% 
        filter(Rsq_CTL_ReHo > 0.1 | Rsq_CTL_fALFF > 0.1)

df <- data.frame(Gene = c(df$Gene,df$Gene,df2$Gene), Direction = c(df$Direction,rep("AllGene",nrow(df)),rep("TopRsq",nrow(df2))))
write.table(df,"integrative_results/ASD_fMRI_Genes.txt",sep="\t",quote=F,row.names=F)
openxlsx::write.xlsx(df, file = "integrative_results/ASD_fMRI_Genes.xlsx", colNames = TRUE, borders = "columns")

###################
## Scatter plot highlighting top 5 genes (pos/neg correlated)
top_labelled <- tbl_df(DiffCor_Sign_Final) %>% 
                  group_by(Direction) %>% 
                  top_n(n = 5, wt = abs(Rho_CTL_fALFF))

pdf("integrative_visualizations/Scatter_Top_Genes.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(DiffCor_Sign_Final, 
            x = "Rho_CTL_ReHo", 
            y = "Rho_CTL_fALFF",
            color = "Direction",
            palette=c("gray60","black"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("ReHo's rho")+ 
      ylab("fALFF's rho")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                      mapping = aes(label = Gene), 
                      size = 5,
                      box.padding = unit(0.4, "lines"),
                      point.padding = unit(0.4, "lines"))+
      theme(legend.position="none")+
      ylim(-0.4,0.4)+
      xlim(-0.4,0.4)
dev.off()

## Boxplot Rsquared
tmpAll <- DiffCor_Sign_Final %>%
          select(Rsq_CTL_fALFF,Rsq_ASD_fALFF) %>%
          melt()

a <- ggboxplot(tmpAll, 
            x = "variable", 
            y = "value",
            color = "variable",
            add = "jitter", 
            line.color = "ligthgrey",
            notch = TRUE, 
            add.params = list(size = 1,alpha=0.3),
            palette = c("black", "red"))+
        stat_compare_means(paired = TRUE,label.x = 1.5,label = "p.signif")+
        xlab("")+ 
        ylab("R^2")+
        theme(legend.position="none") +
        ggtitle("fALFF")+
        rotate_x_text(angle = 45) +
        scale_x_discrete(labels= c("CTL","ASD"))


tmpAll <- DiffCor_Sign_Final %>%
          select(Rsq_CTL_ReHo,Rsq_ASD_ReHo) %>%
          melt()

b <- ggboxplot(tmpAll, 
            x = "variable", 
            y = "value",
            color = "variable",
            add = "jitter", 
            notch = TRUE,
            line.color = "ligthgrey", 
            add.params = list(size = 1,alpha=0.3),
            palette = c("black", "red"))+
        stat_compare_means(paired = TRUE,label.x = 1.5,label = "p.signif")+
        xlab("")+ 
        ylab("R^2")+
        theme(legend.position="none") +
        ggtitle("ReHo")+
        rotate_x_text(angle = 45) +
        scale_x_discrete(labels= c("CTL","ASD"))

plot2by2 <- cowplot::plot_grid(a,b,labels=c("A", "B"), ncol = 2)
cowplot::save_plot("integrative_visualizations/Rsquared_Boxplot.pdf", plot2by2, ncol = 2,base_height=3,base_width=3)


## Compare correlations between ASD and CTL ##
tmpR <- DiffCor_Sign_Final %>% 
            select(Gene,Rho_CTL_ReHo,Rho_ASD_ReHo) %>%
            melt()

A <- ggplot(tmpR, aes(x=value, y=variable, color=variable, point_color=variable,  point_shape = variable, point_size = abs(value),fill=variable)) +
  geom_density_ridges(
          jittered_points=TRUE, 
          alpha = 0.25, 
          point_alpha = 0.5, 
          point_size = 1, 
          position = "raincloud",
          scale=2) +
      scale_fill_manual(values = c("black", "red"), labels = c("Rho_CTL_ReHo", "Rho_ASD_ReHo")) +
      scale_color_manual(values = c("black", "red"), guide = "none") +
      scale_discrete_manual("point_color", values = c("black", "red"), guide = "none") +
      scale_y_discrete(labels= c("CTL","ASD")) +
      guides(fill = guide_legend(
        override.aes = list(
        fill = c("black", "red"),
        color = NA, point_color = NA))
      ) +
      ggtitle("ReHo") +
      theme_ridges(center = TRUE)+
      theme(legend.position="none")+
      xlim(-0.5,0.5)+
      xlab("ReHo's rho")+ 
      ylab("")


tmpR <- DiffCor_Sign_Final %>% 
            select(Gene,Rho_CTL_fALFF,Rho_ASD_fALFF) %>%
            melt() 

B <- ggplot(tmpR, aes(x=value, y=variable, color=variable, point_color=variable,  point_shape = variable, point_size = abs(value),fill=variable)) +
  geom_density_ridges(
          jittered_points=TRUE, 
          alpha = 0.25, 
          point_alpha = 0.5,
          point_size = 1, 
          position = "raincloud",
          scale=2) +
      scale_fill_manual(values = c("black", "red"), labels = c("Rho_CTL_fALFF", "Rho_ASD_fALFF")) +
      scale_color_manual(values = c("black", "red"), guide = "none") +
      scale_discrete_manual("point_color", values = c("black", "red"), guide = "none") +
      scale_y_discrete(labels= c("CTL","ASD")) +
      guides(fill = guide_legend(
        override.aes = list(
        fill = c("black", "red"),
        color = NA, point_color = NA))
        ) +
      ggtitle("fALFF") +
      theme_ridges(center = TRUE)+
      theme(legend.position="none")+
      xlim(-0.5,0.5)+
      xlab("fALFF's rho")+ 
      ylab("")

plot2by2 <- cowplot::plot_grid(A,B,labels=c("A", "B"), ncol = 2)
cowplot::save_plot("integrative_visualizations/Density_Plot_Correlation.pdf", plot2by2, ncol = 2,base_height=3,base_width=3)

# Scatter Rsquared.
# Highlight the genes with more than > 0.1. 
tmp <- DiffCor_Sign_Final %>%
              mutate(Cool = case_when(Rsq_CTL_fALFF > 0.1 & Rsq_CTL_ReHo > 0.1 ~ "Both_High",Rsq_CTL_ReHo > 0.1 ~ "ReHo_High", Rsq_CTL_fALFF > 0.1 ~ "fALFF_High"))
tmp[is.na(tmp)] <- "Others"

top_labelled <- tbl_df(tmp) %>% 
                  filter(Cool == "Both_High") %>% 
                  top_n(n = 6, wt = abs((Rsq_CTL_ReHo+Rsq_CTL_fALFF)/2))

pdf("integrative_visualizations/Scatter_Rsquared.pdf",width=5,height=5,useDingbats=FALSE)
ggscatter(tmp, 
            x = "Rsq_CTL_ReHo", 
            y = "Rsq_CTL_fALFF",
            color = "Cool",
            palette = c("gold","green","gray60","magenta"),
            size = 2,
            alpha=0.5,
            shape=19)+
      xlab("R^2 ReHo")+ 
      ylab("R^2 fALFF")+
      geom_vline(xintercept = 0.1, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 0.1, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      theme(legend.position="none")+
      ylim(0,0.2)+
      xlim(0,0.2) +
      geom_text_repel(data = top_labelled, 
                      mapping = aes(label = Gene), 
                      size = 5,
                      box.padding = unit(0.4, "lines"),
                      point.padding = unit(0.4, "lines"))
dev.off()

## Vulcano plot ##
dge <- read.table(here("dge","DGE_ASDvsCTL.txt"),sep="\t",header=T) # load expression data
dge <- dge[rownames(dge) %in% DiffCor_Sign_Final$Gene,]

tmp <- dge %>%
    rownames_to_column("Names") %>%                    
    mutate(FDR = -log10(adj.P.Val),
           Threshold = case_when(adj.P.Val < 0.05 ~ "DGE", adj.P.Val > 0.05 ~ "NotDGE"),
           Direction = case_when(logFC > 0 ~ "UP",logFC < 0 ~ "DOWN"))


top_labelled <- tbl_df(tmp[tmp$adj.P.Val < 0.05,]) %>% group_by(Direction) %>% top_n(n = 5, wt = FDR)

pdf("integrative_visualizations/Vulcano_Plot_CoolGenes_Expression.pdf",width=5,height=5,useDingbats=FALSE)
ggscatter(tmp, 
      x = "logFC", 
      y = "FDR", 
      color = "Threshold",
      shape = 19,
      size=1,
      palette = c("cyan4", "grey60"))+ 
    geom_hline(yintercept = 1.3, colour = "red",linetype="dotted",size=1,alpha=0.5) +
    geom_vline(xintercept = 0, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
    geom_text_repel(data = top_labelled, mapping = aes(label = Names), size = 5,box.padding = unit(0.4, "lines"),point.padding = unit(0.4, "lines"))+
    theme(legend.position="none")+
    xlim(-1,1)+
    xlab("log2(FC)")+ 
    ylab("-log10(FDR)")
dev.off()

# Upset Plot with previoys study gene set
ctl <- data.frame(Gene = DiffCor_Sign_Final$Gene, Class = rep("DiffCor",nrow(DiffCor_Sign_Final)))
Neuron <- read.table("rawdata/fALFF_genes_Pearson.txt",header=T)
Neuron <- data.frame(Gene=Neuron$fALFF_genes, Class=rep("Neuron",nrow(Neuron)))
Rich <- read.table("rawdata/Richiardi_science.txt",header=T)

tmpUpset <- rbind(ctl,Neuron,Rich)
pdf("integrative_visualizations/Intersection_UpSet.pdf",width=4,height=3,useDingbats=FALSE)
l <- split(as.character(tmpUpset$Gene),tmpUpset$Class)
Class <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Class, ToTGene))
names(metadata) <- c("Class", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))
upset(fromList(l),,nsets = 6, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), list(type = "matrix_rows", 
    column = "sets", colors = c(CTL = "steelblue", Neuron = "purple", Rich = "green"), 
    alpha = 0.5))))
dev.off()

ctl <- data.frame(Gene = CTL_Sign_Final$Gene, Class = rep("CTL",nrow(CTL_Sign_Final)))
Neuron <- read.table("rawdata/fALFF_genes_Pearson.txt",header=T)
Neuron <- data.frame(Gene=Neuron$fALFF_genes, Class=rep("Neuron",nrow(Neuron)))
Rich <- read.table("rawdata/Richiardi_science.txt",header=T)

tmpUpset <- rbind(ctl,Neuron,Rich)
pdf("integrative_visualizations/Intersection_CTL_UpSet.pdf",width=4,height=3,useDingbats=FALSE)
l <- split(as.character(tmpUpset$Gene),tmpUpset$Class)
Class <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Class, ToTGene))
names(metadata) <- c("Class", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))
upset(fromList(l),,nsets = 6, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), list(type = "matrix_rows", 
    column = "sets", colors = c(CTL = "steelblue", Neuron = "purple", Rich = "green"), 
    alpha = 0.5))))
dev.off()

# Coexpression
tmpExpAsd <- p_asd[rownames(p_asd)%in%DiffCor_Sign_Final$Gene,]
tmpExpCtl <- p_ctl[rownames(p_ctl)%in%DiffCor_Sign_Final$Gene,]

coexpASD <- cor(t(tmpExpAsd),method="spearman")
coexpCTL <- cor(t(tmpExpCtl),method="spearman")

pdf("integrative_visualizations/Coexpression_matrix_ASD.pdf",width=4,height=3,useDingbats=FALSE)
pheatmap::pheatmap(
  mat               = coexpASD,
  color             = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
)
dev.off()

pdf("integrative_visualizations/Coexpression_matrix_CTL.pdf",width=4,height=3,useDingbats=FALSE)
pheatmap::pheatmap(
  mat               = coexpCTL,
  color             = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
)
dev.off()

coexpASD <- coexpASD %>%
              as.data.frame() %>%
               rownames_to_column('Gene') %>%
               melt()
coexpASD$Class <- rep("ASD",nrow(coexpASD))

coexpCTL <- coexpCTL %>%
             as.data.frame() %>%
               rownames_to_column('Gene') %>%
               melt()
coexpCTL$Class <- rep("CTL",nrow(coexpCTL))

histCoexp <- rbind(coexpASD,coexpCTL)
histCoexp$abs <- abs(histCoexp$value)
histCoexp$abs[histCoexp$abs == 1] <- NA
histCoexp <- na.omit(histCoexp)

pdf("integrative_visualizations/Coexpression_matrix_ASD_CTL_Histogram.pdf",width=4,height=2,useDingbats=FALSE)
gghistogram(histCoexp, 
  x = "abs",
  add = "mean", 
  fill = "Class",
  palette = c("red", "black"),
  alpha=0.5)+
theme_classic()+
theme(legend.position="right")+
labs(title="",x="abs(rho)", y = "Count")
dev.off()

pdf("integrative_visualizations/Coexpression_matrix_ASD_CTL_Histogram_2.pdf",width=4,height=2,useDingbats=FALSE)
gghistogram(histCoexp, 
  x = "value",
  fill = "Class",
  palette = c("red", "black"),
  alpha=0.5)+
theme_classic()+
theme(legend.position="right")+
labs(title="",x="rho", y = "Count")
dev.off()

# DGCA on TopRsq
genes <- read.table("integrative_results/ASD_fMRI_Genes.txt",header=T)
genes <- genes %>% filter(Direction == "TopRsq")

tmpExpAsd <- p_asd[rownames(p_asd)%in%genes$Gene,]
tmpExpCtl <- p_ctl[rownames(p_ctl)%in%genes$Gene,]
exp <- cbind(tmpExpCtl,tmpExpAsd)

tmpdemo <- rbind(demoCTL,demoASD)
tmpdemo <- tmpdemo %>%
                        mutate(CTL= case_when(Diagnosis == "CTL" ~ 1, Diagnosis == "ASD"  ~ 0),
                               ASD= case_when(Diagnosis == "CTL" ~ 0, Diagnosis == "ASD"  ~ 1))
design_mat=as.matrix(tmpdemo[,c("ASD","CTL")])

pdf("integrative_visualizations/Coexpression_matrix_ASD_CTL_DGCA.pdf",width=10,height=8,useDingbats=FALSE)
ddcor_res = ddcorAll(inputMat = exp, design = design_mat,classify = TRUE,compare = c("ASD", "CTL"),corrType = "spearman",adjust = "bonferroni", nPerms = 100, nPairs = "all",heatmapPlot = TRUE)
dev.off()

openxlsx::write.xlsx(ddcor_res, file = "integrative_results/DGCA_Output.xlsx", colNames = TRUE, borders = "columns",sheetName="DGCA")

# fold change heatmap

load(here("dge", "DGE_ASDvsCTL_RegByReg.RData"))
genes <- read.table(here("integrative_results","ASD_fMRI_Genes.txt"),header=T,sep="\t")

genes <- genes %>%
          filter(Direction == "AllGene")

new_names <- names(brDGE)

pdf("integrative_visualizations/FoldChange_Heatmap.pdf", width=3, height=4)
map(brDGE, ~ (.x %>% rownames_to_column("Gene") %>% select(Gene,logFC))) %>% 
      map2(new_names, ~setnames(.x, 'logFC', .y)) %>%
      purrr::reduce(left_join, by = "Gene") %>% 
      filter(Gene %in% genes$Gene) %>% 
      column_to_rownames("Gene") %>% 
      pheatmap::pheatmap(
          color         = colorRampPalette(c("navy", "white", "firebrick3"))(50),
          border_color  = NA,
          show_colnames = TRUE,
          show_rownames = FALSE,
          scale="row"
        )
dev.off()
      
# Venn diagram Development TopRsq
dev <- read.table("integrative_results/ASD_fMRI_DevClustered.txt",header=T)
all <- read.table("integrative_results/ASD_fMRI_Genes.txt",header=T)
all <- all %>%
    filter(Direction == "TopRsq") %>%
    droplevels() %>%
    rename(Clusters = Direction)

tmp <- rbind(dev,all)
tmpListed <- split(tmp,tmp$Clusters)

# salmon = ff8c69 
# olivedrab = 6b8e23
# thistle = d8bfd8

ven <- venn.diagram(x= list(
Adult = tmpListed[[1]]$Gene,
EarlyDev = tmpListed[[2]]$Gene,
Stable = tmpListed[[3]]$Gene,
TopRsq = tmpListed[[4]]$Gene),
filename = NULL,
fill = c("salmon","olivedrab","thistle","gold"),
alpha = 0.50, 
cex = 1, 
fontfamily = "serif", 
cat.cex = 1.5,
cat.fontfamily = "serif",
cat.dist = 0.2,
margin = 0.2,
lty = "dashed")

pdf(file="integrative_visualizations/Venn_Diagram_Dev_TopRsq.pdf",width=4,height=4)
grid.draw(ven)
dev.off()

