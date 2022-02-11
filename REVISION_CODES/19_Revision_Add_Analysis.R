rm(list=ls())
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pls))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(effsize))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rio))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(STRINGdb))
suppressPackageStartupMessages(library(igraph))
source("UTILS/Utils.R")

# Create Output directories
folder_names <- c("revision")
sapply(folder_names, dir.create)


load("integrative_results/DiffCor_Significant.RData")
load("integrative_results/DiffCor_Tables.RData")

df_sign <- DiffCor_Sign_Final %>%
      select(Gene,Rsq_CTL_fALFF, Rsq_CTL_ReHo,DiffCor_ReHo_Z,DiffCor_fALFF_Z,DiffCor_Comb_P,Direction) %>%
      mutate(Mean_Z = (DiffCor_ReHo_Z+DiffCor_fALFF_Z)/2, Mean_Rsq = (Rsq_CTL_ReHo+Rsq_CTL_fALFF)/2) %>%
      mutate(Rank_Z = rank(Mean_Z),Rank_Rsq = rank(Mean_Rsq))

min <- min(df_sign$Mean_Z)


df <- DiffCor_Combined %>%
      select(Gene,Rsq_CTL_fALFF, Rsq_CTL_ReHo,DiffCor_ReHo_Z,DiffCor_fALFF_Z,DiffCor_Comb_P) %>%
      mutate(Mean_Z = (DiffCor_ReHo_Z+DiffCor_fALFF_Z)/2, Mean_Rsq = (Rsq_CTL_ReHo+Rsq_CTL_fALFF)/2) %>%
      mutate(Rank_Z = rank(Mean_Z),Rank_Rsq = rank(Mean_Rsq))


top_labelled <- as_tibble(df) %>% 
                  top_n(n = 5, wt = abs(Rank_Z))
                  #filter(Gene %in% c("SCN1B","PVALB","SYT2","FILIP1"))

pdf("revision/Effect_Size_Plot.pdf",width=7,height=4)
key.gns <- c("GABRQ", "SCN1B", "PVALB","FILIP1")
ggdensity(df, x = "Mean_Z", y = "..count..",
            xlab = "Differential Correlation Z",
            ylab = "Number of genes",
            fill = "lightgray", color = "black",
            label = "Gene", label.select = key.gns, repel = TRUE,
            font.label = list(color= "Mean_Z"),
            xticks.by = 0.5, # Break x ticks by 20
            gradient.cols = c("blue", "red"),
            legend = c(0.7, 0.6),                                 
            legend.title = ""       # Hide legend title
            ) +
  geom_vline(xintercept=min, colour="grey") +
  geom_text(aes(x=min, label="\nSignificant", y=9000), colour="blue", angle=90, text=element_text(size=12))
dev.off()


# Pvalb cor vs PCA
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

# Loading imaging data
load(here("rawdata","Imaging_Normalized.RData"))

# Reorder the expression data
p_ctl1 <- p_ctl[,match(do.call(rbind,split(demoCTL,demoCTL$Region))$ID,colnames(p_ctl))]
p_asd1 <- p_asd[,match(do.call(rbind,split(demoASD,demoASD$Region))$ID,colnames(p_asd))]


#############################
### Pan Cortical Analysis ###
#############################
B=200  ## select number of times for leave-multiple-out method
tmpCTL_PanCort <- vector("list", length = B)
names <- unique(do.call(rbind,split(demoCTL,demoCTL$Region))$Region)
for (i in 1:B)
  { 
  tmpCTL_PanCort[[i]] <- Imaging %>%
                            group_by(Region) %>%  
                            nest() %>%           
                            arrange(Region, names) %>% 
                            add_column(n = sapply(split(demoCTL,demoCTL$Region), nrow)) %>% 
                            mutate(samp = map2(data, n, sample_n)) %>% 
                            select(Region, samp) %>%
                            unnest() %>%
                            as.data.frame()
}

B=200  ## select number of times for leave-multiple-out method
tmpASD_PanCort <- vector("list", length = B)
for (i in 1:B)
  {
  tmpASD_PanCort[[i]] <- Imaging %>%
                            group_by(Region) %>% 
                            nest() %>%
                            arrange(Region, names) %>%            
                            add_column(n = sapply(split(demoASD,demoASD$Region), nrow)) %>% 
                            mutate(samp = map2(data, n, sample_n)) %>% 
                            select(Region, samp) %>% 
                            unnest() %>%
                            as.data.frame()
}


# Analysis PCA vs fALFF1 or ReHo1
pca.ctl<-prcomp(t(p_ctl1))

PCi<-data.frame(pca.ctl$x[,1:2],fALFF1=tmpCTL_PanCort[[1]]$fALFF1,ReHo = tmpCTL_PanCort[[1]]$ReHo1, PVALB = p_ctl1["PVALB",])


a <- ggscatter(PCi, x = "PC1", y = "fALFF1",
          add = "reg.line",                         # Add regression line
          color = "darkgreen",                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )+
  stat_cor()           # Add correlation coefficient


b <- ggscatter(PCi, x = "PC2", y = "fALFF1",
          add = "reg.line",                         # Add regression line
          color = "darkgreen",                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )+
  stat_cor() 


c <- ggscatter(PCi, x = "PC1", y = "ReHo",
          add = "reg.line",                         # Add regression line
          color = "magenta",                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )+
  stat_cor()           # Add correlation coefficient


d <- ggscatter(PCi, x = "PC2", y = "ReHo",
          add = "reg.line",                         # Add regression line
          color = "magenta",                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )+
  stat_cor() 


e <- ggscatter(PCi, x = "PVALB", y = "fALFF1",
          add = "reg.line",                         # Add regression line
          color = "blue",                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )+
  stat_cor()           # Add correlation coefficient


f <- ggscatter(PCi, x = "PVALB", y = "ReHo",
          add = "reg.line",                         # Add regression line
          color = "blue",                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )+
  stat_cor() 




plot <- plot_grid(e,f,a,b,c,d, labels=c("A", "B","C","D","E","F"), ncol = 2,nrow=3)
save_plot("revision/fMRIvsPCs.pdf", plot, ncol = 2,base_height=8,base_width=4)


# Scatter  

tmp_a<-data.frame(fALFF1=tmpCTL_PanCort[[1]]$fALFF1,
                  ReHo = tmpCTL_PanCort[[1]]$ReHo1, 
                  PVALB = p_ctl1["PVALB",], 
                  SCN1B = p_ctl1["SCN1B",],
                  FILIP1 = p_ctl1["FILIP1",],
                  GABRQ = p_ctl1["GABRQ",],
                  Class = "CTL")

tmp_b<-data.frame(fALFF1=tmpASD_PanCort[[1]]$fALFF1,
                  ReHo = tmpASD_PanCort[[1]]$ReHo1, 
                  PVALB = p_asd1["PVALB",], 
                  SCN1B = p_asd1["SCN1B",],
                  FILIP1 = p_asd1["FILIP1",],
                  GABRQ = p_asd1["GABRQ",],
                  Class = "ASD")

df_scatter <- rbind(tmp_a, tmp_b)

a <-  ggscatter(df_scatter, x = "GABRQ", y = "fALFF1",
          add = "reg.line",                         # Add regression line
          color = "Class",
          palette = c("red","black"),                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )

b <-  ggscatter(df_scatter, x = "FILIP1", y = "fALFF1",
          add = "reg.line",                         # Add regression line
          color = "Class",
          palette = c("red","black"),                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )


c <-  ggscatter(df_scatter, x = "PVALB", y = "fALFF1",
          add = "reg.line",                         # Add regression line
          color = "Class",
          palette = c("red","black"),                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )


d <-  ggscatter(df_scatter, x = "SCN1B", y = "fALFF1",
          add = "reg.line",                         # Add regression line
          color = "Class",
          palette = c("red","black"),                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )

plot <- plot_grid(a,b,c,d, labels=c("GABRQ", "FILIP1","PVALB","SCN1B"), ncol = 2,nrow=2)
save_plot("revision/Example_DiffCor_Genes.pdf", plot, ncol = 2,base_height=8,base_width=4)
