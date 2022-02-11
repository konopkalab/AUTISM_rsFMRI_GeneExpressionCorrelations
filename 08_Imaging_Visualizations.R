rm(list=ls())
suppressPackageStartupMessages({
library(ggplot2)
library(ggjoy)
library(ggrepel)
library(ggpubr)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(cowplot)
library(GGally)
library(ggforce)
library(VennDiagram)
library(ggstatsplot)
library(broom)
library(tidyverse)
library(here)
library(effsize)
})

# Create Directory
folder_names <- c("imaging_visualizations")
sapply(folder_names, dir.create)

# PCA plot
fALFF1 <- read.table(here("rawdata","ABIDEIandIIChosenScansBrodmannValues_equalized_global_fALFF.csv"),sep=",",header=T)
colnames(fALFF1) <- gsub("MeanRegion","BA",colnames(fALFF1))

fALFF1 <- fALFF1 %>% 
      mutate(Diagnosis = case_when(DX == 1 ~ "ASD", DX == 2 ~ "CTL"),
             Sex = case_when(SEX == 2 ~ "F", SEX == 1 ~ "M")) %>%
      select(-DX, -SEX) %>%
      select(Sex, everything()) %>%
      select(Diagnosis, everything()) %>%
      select(ID, everything()) %>%
      as.data.frame() 

PCA <- prcomp(fALFF1[,grep("BA",colnames(fALFF1))])

pdf("imaging_visualizations/PCA_fALFF_ScreePlot.pdf",width=6,height=4)
factoextra::fviz_eig(PCA,ncp = 25)
dev.off()

PCi<-data.frame(PCA$x,Diagnosis=fALFF1$Diagnosis)
eig <- (PCA$sdev)^2
variance <- eig*100/sum(eig)

pdf("imaging_visualizations/PCA_fALFF.pdf",width=5,height=4)
ggscatter(PCi, x = "PC1", y = "PC2",color = "Diagnosis",palette=c("red","black"), shape = 21, size = 1)+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic()
dev.off()

ReHo1 <- read.table(here("rawdata","ABIDEIandIIChosenScansBrodmannValues_equalized_global_ReHo.csv"),sep=",",header=T)
colnames(ReHo1) <- gsub("MeanRegion","BA",colnames(ReHo1))

ReHo1 <- ReHo1 %>% 
      mutate(Diagnosis = case_when(DX == 1 ~ "ASD", DX == 2 ~ "CTL"),
             Sex = case_when(SEX == 2 ~ "F", SEX == 1 ~ "M")) %>%
      select(-DX, -SEX) %>%
      select(Sex, everything()) %>%
      select(Diagnosis, everything()) %>%
      select(ID, everything()) %>%
      as.data.frame() 

PCA <- prcomp(ReHo1[,grep("BA",colnames(ReHo1))])

pdf("imaging_visualizations/PCA_ReHo_ScreePlot.pdf",width=6,height=4)
factoextra::fviz_eig(PCA,ncp = 25)
dev.off()

PCi<-data.frame(PCA$x,Diagnosis=ReHo1$Diagnosis)
eig <- (PCA$sdev)^2
variance <- eig*100/sum(eig)

pdf("imaging_visualizations/PCA_ReHo.pdf",width=5,height=4)
ggscatter(PCi, x = "PC1", y = "PC2",color = "Diagnosis",palette=c("red","black"), shape = 21, size = 1)+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic()
dev.off()


# Viz of the ABIDE values with ASD
load(here("rawdata","Imaging_with_ASD_Normalized.RData"))
#Imaging_ASD$Region <- factor(Imaging_ASD$Region,levels = c("BA9","BA44_45","BA24","BA38","BA20_37","BA4_6","BA3_1_2_5","BA39_40","BA41_42_22","BA7","BA17"))

# Cohen's d calcualtion
tmpCast <- split(Imaging_ASD,Imaging_ASD$Region)
cohenF <- tmpCast %>% 
        purrr::map(~ cohen.d(.x$fALFF1,.x$Diagnosis)$estimate) %>%
        purrr::map_df(tidy) %>%
        mutate(Region = names(tmpCast)) %>%
        as.data.frame() %>%
        dplyr::rename(fALFF = x)
        
openxlsx::write.xlsx(cohenF, file = "imaging_visualizations/CohenF_fALFF_Stats.xlsx", colNames = TRUE, borders = "columns",overwrite=TRUE)

cohenR <- tmpCast %>% 
        map(~ cohen.d(.x$ReHo1,.x$Diagnosis)$estimate) %>%
        map_df(tidy) %>%
        mutate(Region = names(tmpCast)) %>%
        dplyr::rename(ReHo = x) %>%
        as.data.frame()
openxlsx::write.xlsx(cohenR, file = "imaging_visualizations/CohenF_ReHo_Stats.xlsx", colNames = TRUE, borders = "columns",overwrite=TRUE)

pdf(file="imaging_visualizations/Imaging_CoehnsD.pdf",width=6,height=3)
full_join(cohenF,cohenR) %>%
    melt() %>%
    ggplot(aes(Region, variable, fill= value)) + 
        geom_tile()+
        scale_fill_gradientn(colours=c('red','white','blue'), limits=c(-1,1)) +
        geom_text(aes(Region, variable, label = round(value,digits=2)), color = "black", size = 3)+
          theme_classic()+
          rotate_x_text(angle = 45)+
          xlab("")+
          ylab("")+ 
          labs(fill = "Cohen's d")
dev.off()

# Boxplotting
pdf(file="imaging_visualizations/Imaging_boxplot_fALFF_Simplified.pdf",width=8,height=4)
tmpMelted <- Imaging_ASD %>% melt() %>% filter(variable == "fALFF1") %>% arrange(Region)
ggboxplot(tmpMelted, 
          x = "Region", 
          y = "value",
          color = "Diagnosis", 
          palette = c("red","black"), 
          short.panel.labs = FALSE)+
          #facet_wrap(~Region,ncol=4,nrow=3,scales="free")+
  stat_compare_means(aes(group = Diagnosis),label = "p.signif",hide.ns = TRUE)+
theme_classic()+
theme(legend.position="none")+
rotate_x_text(angle = 45)+
xlab("")+
ylab("fALFF")
dev.off()

Stats <- compare_means(value ~ Diagnosis, tmpMelted, group.by = "Region",method="t.test") %>% as.data.frame()
openxlsx::write.xlsx(Stats, file = "imaging_visualizations/Imaging_boxplot_fALFF_Stats.xlsx", colNames = TRUE, borders = "columns",overwrite=TRUE)

pdf(file="imaging_visualizations/Imaging_boxplot_ReHo_Simplified.pdf",width=8,height=4)
tmpMelted <- Imaging_ASD %>% melt() %>% filter(variable == "ReHo1") %>% arrange(Region)
ggboxplot(tmpMelted, 
          x = "Region", 
          y = "value",
          color = "Diagnosis", 
          palette = c("red","black"), 
          short.panel.labs = FALSE)+
          #facet_wrap(~Region,ncol=4,nrow=3,scales="free")+
  stat_compare_means(aes(group = Diagnosis),label = "p.signif", hide.ns = TRUE)+
theme_classic()+
theme(legend.position="none")+
rotate_x_text(angle = 45)+
xlab("")+
ylab("ReHo")
dev.off()

Stats <- compare_means(value ~ Diagnosis, tmpMelted, group.by = "Region",method="t.test") %>% as.data.frame()
openxlsx::write.xlsx(Stats, file = "imaging_visualizations/Imaging_boxplot_ReHo_Stats.xlsx", colNames = TRUE, borders = "columns",overwrite=TRUE)


# Working on CTL only
load(here("rawdata","Imaging_Normalized.RData"))

# Set up colors for each region
paired <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00","#CAB2D6","#6A3D9A","#FF99F8")

pdf(file="imaging_visualizations/Imaging_ScatterCor.pdf",width=7,height=6)
ggscatter(Imaging, 
          x = "fALFF1", 
          y = "ReHo1",
          add = "reg.line",                         # Add regression line
          color = "Region",
          palette=paired,            # Color by groups "cyl"
          fullrange = FALSE,                         # Extending the regression line
          rug = FALSE,
          ellipse = TRUE,size = 1                                 # Add marginal rug
          )+
  stat_cor(aes(color = Region)) +
  theme(legend.position="right")
dev.off()

svg("imaging_visualizations/Imaging_ScatterCor.svg")
ggscatter(Imaging, 
          x = "fALFF1", 
          y = "ReHo1",
          add = "reg.line",                         # Add regression line
          color = "Region",
          palette=paired,            # Color by groups "cyl"
          fullrange = FALSE,                         # Extending the regression line
          rug = FALSE,
          ellipse = TRUE,size = 1                                 # Add marginal rug
          )+
  stat_cor(aes(color = Region)) +
  theme(legend.position="right")
dev.off()


pdf(file="imaging_visualizations/Imaging_ScatterCor_PanCort.pdf",width=6,height=6)
ggscatter(Imaging, x = "fALFF1", y = "ReHo1",
   color = "Region", size = 1, 
   add = "reg.line", palette=paired, 
   add.params = list(color = "red", fill = "lightgray"), 
   conf.int = TRUE, 
   cor.coef = TRUE, 
   cor.coeff.args = list(method = "pearson", label.sep = "\n"),
   ellipse = TRUE)+
  theme(legend.position="right")
dev.off()

# PCA plots
pdf(file="imaging_visualizations/PCA_fALFF.pdf",width=5,height=5)
falff <- dcast(Imaging,ID ~ Region,value.var = "fALFF1")
falff$ID <- NULL
PC<-prcomp(t(falff))
eig <- (PC$sdev)^2
variance <- eig*100/sum(eig)
PCi<-data.frame(PC$x[,1:3],Class = rownames(PC$x))
ggscatter(PCi, x = "PC1", y = "PC2",color = "Class",
palette = paired,
ellipse = FALSE, mean.point = FALSE,star.plot = FALSE,rug = FALSE,label="Class")+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )")) +
theme(legend.position="none")+
border(color = "green")
dev.off()


pdf(file="imaging_visualizations/PCA_ReHo.pdf",width=5,height=5)
reho <- dcast(Imaging,ID ~ Region,value.var = "ReHo1")
reho$ID <- NULL
PC<-prcomp(t(reho))
eig <- (PC$sdev)^2
variance <- eig*100/sum(eig)
PCi<-data.frame(PC$x[,1:3],Class = rownames(PC$x))
ggscatter(PCi, x = "PC1", y = "PC2",color = "Class",
palette = paired,
ellipse = FALSE, mean.point = FALSE,star.plot = FALSE,rug = FALSE,label="Class")+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )")) +
theme(legend.position="none")+
border(color = "magenta")
dev.off()

# crosscorr plots
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}


pdf(file="imaging_visualizations/CrossCorr_fALFF.pdf",width=8,height=8)
ggpairs(falff, lower = list(continuous = my_fn))
dev.off()

pdf(file="imaging_visualizations/CrossCorr_ReHo.pdf",width=8,height=8)
ggpairs(reho, lower = list(continuous = my_fn))
dev.off()

pdf(file="imaging_visualizations/CrossCorr_fALFF_2.pdf",width=5,height=5)
ggcorr(falff, palette = "RdBu", label = TRUE)
dev.off()

pdf(file="imaging_visualizations/CrossCorr_ReHo_2.pdf",width=5,height=5)
ggcorr(reho, palette = "RdBu", label = TRUE)
dev.off()


########################### 
### Distribution values ###
###########################
a <- ggboxplot(Imaging, 
          x = "Region", 
          y = "fALFF1",
          color = "Region", 
          palette = paired, 
          short.panel.labs = FALSE)+
    theme_classic()+
    theme(legend.position="none")+
    rotate_x_text(angle = 45)+
    xlab("")+
    ylab("fALFF")

b <- ggboxplot(Imaging, 
          x = "Region", 
          y = "ReHo1",
          color = "Region", 
          palette = paired, 
          short.panel.labs = FALSE)+
    theme_classic()+
    theme(legend.position="none")+
    rotate_x_text(angle = 45)+
    xlab("")+
    ylab("ReHo")

plot <- plot_grid(a,b,labels=c("A", "B"), ncol = 1,nrow=2)
save_plot("imaging_visualizations/Distribution_Imaging_Values.pdf", plot, ncol = 1,base_height=6,base_width=5)

sessionInfo()


