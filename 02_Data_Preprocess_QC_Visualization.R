# Covariate plot
rm(list=ls())
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggjoy))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))

# load data you need for visualization
load(here("rawdata","Expression_Data.RData"))

################## 
###   VarExp   ###
##################
#var <- VarExp(exp,demo,2,FALSE)
#pdf("PLOTS/Variance_Explained.pdf",width=6,height=6)
#plotVarExp(var,"Variance Explained")
#dev.off()

################## 
###  PCA plot  ###
##################
PCA <- prcomp(t(exp))

pdf("rawdata/PCA_ScreePlot.pdf",width=6,height=4)
factoextra::fviz_eig(PCA,ncp = 25)
dev.off()

PCi<-data.frame(PCA$x,Diagnosis=demo$Diagnosis,Region = demo$Region, Labels = rownames(demo),Age = demo$Age, Bank = demo$BrainBank, PMI = demo$PMI, Batch = demo$Batch)
eig <- (PCA$sdev)^2
variance <- eig*100/sum(eig)

pcaExp_Diagnosis <- ggscatter(PCi, x = "PC1", y = "PC2",color = "Diagnosis",palette="jco", shape = 21, size = 3)+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic()

PCA2<- prcomp(t(expAdj))
PCi2<-data.frame(PCA2$x,Diagnosis=demo$Diagnosis,Region = demo$Region, Labels = rownames(demo),Age = demo$Age, Bank = demo$BrainBank, PMI = demo$PMI, Batch = demo$Batch)
eig <- (PCA2$sdev)^2
variance <- eig*100/sum(eig)

pcaAdj_Diagnosis <- ggscatter(PCi2, x = "PC1", y = "PC2",color = "Diagnosis",palette="jco", shape = 21, size = 3)+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic()

plot2by2 <- cowplot::plot_grid(pcaExp_Diagnosis, pcaAdj_Diagnosis,labels=c("A", "B"), ncol = 2)
cowplot::save_plot("rawdata/PCAs.pdf", plot2by2, ncol = 2,base_height=4,base_width=5)

######################### 
### Plot demographics ###
#########################
pdf("rawdata/Demo_Pairs_Demographics.pdf",width=12,height=10)
p <- ggpairs(demo, aes(color = Diagnosis),columns = names(demo)[3:10]) + 
theme_classic()
for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + 
        scale_fill_manual(values=c("steelblue", "orange")) +
        scale_color_manual(values=c("steelblue", "orange"))  
  }
}
print(p)
dev.off()

pdf("rawdata/Demo_Pairs_TechSeq.pdf",width=12,height=10)
p <- ggpairs(demo, aes(color = Diagnosis),columns = names(demo)[11:21]) + 
theme_classic()
for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + 
        scale_fill_manual(values=c("steelblue", "orange")) +
        scale_color_manual(values=c("steelblue", "orange"))  
  }
}
print(p)
dev.off()

######################### 
### Plot demographics ###
#########################
demo$pca1<-PCA$x[,1]
fit1=lm(pca1 ~ .,data=droplevels(demo[c(-2)]))
df <- tidy(fit1)[-1,]
df$log10=-log10(df$p.value)

a <- ggbarplot(df, x = "term", y = "log10",
          fill = "white",           
          color = "blue",            
          x.text.angle = 90,          
          ylab = "-log10(P)",,
          xlab = "Covariates",
          rotate = TRUE,
          ggtheme = theme_classic()
          )+ 
geom_hline(yintercept = 1.3, linetype="dotted", color = "red", size=1)+
ylim(0,16) +
ggtitle("PCA1 vs Cov: Raw")

demo$pca1<-PCA2$x[,1]
fit1=lm(pca1 ~ .,data=droplevels(demo[c(-2)]))
df <- tidy(fit1)[-1,]
df$log10=-log10(df$p.value)

b <- ggbarplot(df, x = "term", y = "log10",
          fill = "white",           
          color = "blue",            
          x.text.angle = 90,          
          ylab = "-log10(P)",
          xlab = "Covariates",
          rotate = TRUE,
          ggtheme = theme_classic()
          )+ 
geom_hline(yintercept = 1.3, linetype="dotted", color = "red", size=1)+
ylim(0,16) +
ggtitle("PCA1 vs Cov: AdjExp")

demo$pca1<-PCA$x[,2]
fit1=lm(pca1 ~ .,data=droplevels(demo[c(-2)]))
df <- tidy(fit1)[-1,]
df$log10=-log10(df$p.value)

c <- ggbarplot(df, x = "term", y = "log10",
          fill = "white",           
          color = "blue",           
          x.text.angle = 90,          
          ylab = "-log10(P)",,
          xlab = "Covariates",
          rotate = TRUE,
          ggtheme = theme_classic()
          )+ 
geom_hline(yintercept = 1.3, linetype="dotted", color = "red", size=1)+
ylim(0,16)+
ggtitle("PCA2 vs Cov: Raw")

demo$pca1<-PCA2$x[,2]
fit1=lm(pca1 ~ .,data=droplevels(demo[c(-2)]))
df <- tidy(fit1)[-1,]
df$log10=-log10(df$p.value)

d <- ggbarplot(df, x = "term", y = "log10",
          fill = "white",           
          color = "blue",           
          x.text.angle = 90,          
          ylab = "-log10(P)",
          xlab = "Covariates",
          rotate = TRUE,
          ggtheme = theme_classic()
          )+ 
geom_hline(yintercept = 1.3, linetype="dotted", color = "red", size=1)+
ylim(0,16)+
ggtitle("PCA2 vs Cov: AdjExp")

plot2by2 <- cowplot::plot_grid(a,b,c,d,labels=c("A", "B","C","D"), ncol = 2)
cowplot::save_plot("rawdata/PCsVsCov.pdf", plot2by2, ncol = 2,base_height=10,base_width=6)

sessionInfo()