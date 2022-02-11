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
library(ggseg)
})


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


# PCA loadings 
pca.ctl<-prcomp(p_ctl1)$rotation

PCi<-data.frame(PC1 = pca.ctl[,1],
				PC2 = pca.ctl[,2],
				fALFF=tmpCTL_PanCort[[1]]$fALFF1,
				ReHo = tmpCTL_PanCort[[1]]$ReHo1,
				PVALB = p_ctl1["PVALB",],
				SCN1B = p_ctl1["SCN1B",],
				FILIP1 = p_ctl1["FILIP1",],
				GABRQ = p_ctl1["GABRQ",])


a <- ggscatter(PCi, x = "PC1", y = "PVALB",
          add = "reg.line",                         # Add regression line
          color = "blue",                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )+
  stat_cor()           # Add correlation coefficient


b <- ggscatter(PCi, x = "PC2", y = "PVALB",
          add = "reg.line",                         # Add regression line
          color = "blue",                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )+
  stat_cor()     


c <- ggscatter(PCi, x = "PC1", y = "SCN1B",
          add = "reg.line",                         # Add regression line
          color = "blue",                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )+
  stat_cor()           # Add correlation coefficient


d <- ggscatter(PCi, x = "PC2", y = "SCN1B",
          add = "reg.line",                         # Add regression line
          color = "blue",                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )+
  stat_cor()

plot <- plot_grid(a,b,c,d,labels=c("a", "","b",""), ncol = 2,nrow=2)
save_plot("revision/PVALB_SCN1BvsPCA.pdf", plot, ncol = 2,base_height=4,base_width=2)


dfBar <- data.frame(
		 Gene = c("PVALB","SCN1B","PVALB","SCN1B","PVALB","SCN1B","PVALB","SCN1B"),
	     PCs = c("PC2","PC2","PC1","PC1","fALFF","fALFF","ReHo","ReHo"), 
		 Cor = c(round(cor(PCi$PC2,PCi$PVALB),2),round(cor(PCi$PC2,PCi$SCN1B),2),
		 		 round(cor(PCi$PC1,PCi$PVALB),2),round(cor(PCi$PC1,PCi$SCN1B),2),
         round(cor(PCi$fALFF,PCi$PVALB),2),round(cor(PCi$fALFF,PCi$SCN1B),2),
         round(cor(PCi$ReHo,PCi$PVALB),2),round(cor(PCi$ReHo,PCi$SCN1B),2)))

pdf("revision/PVALB_SCN1B_Loadings_Correlation.pdf",height=4,width=3)
 ggbarplot(dfBar, "PCs", "Cor",
   fill = "Gene", color = "Gene",
   label = FALSE,
   orientation = "horiz",
   palette = c("#00AFBB", "#E7B800"),
 position = position_dodge(0.9))
dev.off()

# GGseg
df <-data.frame(pca.ctl[,1:2], Region = demoCTL$Region) %>%
  group_by(Region) %>%
  summarise(Loadings = as.numeric(median(PC1)))%>%
  as.data.frame() %>%
  mutate(scaled = scale(Loadings))


someData <- tibble(
  region = c("lateral occipital", "inferior temporal", "caudal anterior cingulate","postcentral","middle temporal","supramarginal",
  			"caudal middle frontal","superior temporal", "pars opercularis","superior parietal","rostral middle frontal"),
  ba =  c("BA17", "BA20_37","BA24","BA3_1_2_5","BA38","BA39_40","BA4_6","BA41_42_22","BA44_45","BA7","BA9"),
  Loadings = df$scaled
)

cortical_pos <- c("left lateral", "left medial")

pdf("revision/Expression_Loadings_GGSeg.pdf",height=5,width=5)
someData %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(cortical_pos),
             aes(fill = Loadings)) +
  scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  theme_void()
dev.off()


df <-data.frame(PVALB = p_ctl1["PVALB",], Region = demoCTL$Region) %>%
  group_by(Region) %>%
  summarise(PVALB_mean = as.numeric(mean(PVALB)))%>%
  as.data.frame() %>%
  mutate(scaled = scale(PVALB_mean))


someData <- tibble(
  region = c("lateral occipital", "inferior temporal", "caudal anterior cingulate","postcentral","middle temporal","supramarginal",
  			"caudal middle frontal","superior temporal", "pars opercularis","superior parietal","rostral middle frontal"),
  ba =  c("BA17", "BA20_37","BA24","BA3_1_2_5","BA38","BA39_40","BA4_6","BA41_42_22","BA44_45","BA7","BA9"),
  PVALB = df$scaled
)

cortical_pos <- c("left lateral", "left medial")

pdf("revision/Expression_PVALB_GGSeg.pdf",height=5,width=5)
someData %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(cortical_pos),
             aes(fill = PVALB)) +
  scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  theme_void()
dev.off()


df <-data.frame(SCN1B = p_ctl1["SCN1B",], Region = demoCTL$Region) %>%
  group_by(Region) %>%
  summarise(SCN1B_mean = as.numeric(mean(SCN1B)))%>%
  as.data.frame() %>%
  mutate(scaled = scale(SCN1B_mean))


someData <- tibble(
  region = c("lateral occipital", "inferior temporal", "caudal anterior cingulate","postcentral","middle temporal","supramarginal",
  			"caudal middle frontal","superior temporal", "pars opercularis","superior parietal","rostral middle frontal"),
  ba =  c("BA17", "BA20_37","BA24","BA3_1_2_5","BA38","BA39_40","BA4_6","BA41_42_22","BA44_45","BA7","BA9"),
  SCN1B = df$scaled
)

cortical_pos <- c("left lateral", "left medial")

pdf("revision/Expression_SCN1B_GGSeg.pdf",height=5,width=5)
someData %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(cortical_pos),
             aes(fill = SCN1B)) +
  scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  theme_void()
dev.off()

# Scatter loadings/Cohens'd 

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
        
cohenR <- tmpCast %>% 
        map(~ cohen.d(.x$ReHo1,.x$Diagnosis)$estimate) %>%
        map_df(tidy) %>%
        mutate(Region = names(tmpCast)) %>%
        dplyr::rename(ReHo = x) %>%
        as.data.frame()

PCscatt<-data.frame(PC1 = pca.ctl[,1],
        PC2 = pca.ctl[,2],
        CohenD_fALFF =cohenF$fALFF,
        CohenD_ReHo =cohenF$ReHo)

df <-data.frame(pca.ctl[,1:2], Region = demoCTL$Region) %>%
  group_by(Region) %>%
  summarise(Loadings = as.numeric(median(PC1)))%>%
  as.data.frame() %>%
  mutate(scaled = scale(Loadings), ,
        CohenD_fALFF =cohenF$fALFF,
        CohenD_ReHo =cohenR$ReHo)

a <- ggscatter(df, x = "Loadings", y = "CohenD_fALFF",
          add = "reg.line",                         # Add regression line
          color = "green",                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )+
  stat_cor()   


b <- ggscatter(df, x = "Loadings", y = "CohenD_ReHo",
          add = "reg.line",                         # Add regression line
          color = "magenta",                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = FALSE                                # Add marginal rug
          )+
  stat_cor()   

