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

pdf(file="revision/Correlation_CoehnsD.pdf",width=4,height=4)
full_join(cohenF,cohenR) %>%
			ggscatter( 
            x = "ReHo", 
            y = "fALFF",
            palette="grey60",
            alpha=0.3,
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = FALSE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = -0.2,lable.y=0.2, label.sep = "\n"))+
      xlab("Cohen's d ReHo")+ 
      ylab("Cohen's d fALFF")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      theme(legend.position="none")+
      ylim(-0.2,0.2)+
      xlim(-0.2,0.2)
dev.off()

# Seg Viz
someData <- tibble(
  region = c("lateral occipital", "inferior temporal", "caudal anterior cingulate","postcentral","middle temporal","supramarginal",
  			"caudal middle frontal","superior temporal", "pars opercularis","superior parietal","rostral middle frontal"),
  ba =  c("BA17", "BA20_37","BA24","BA3_1_2_5","BA38","BA39_40","BA4_6","BA41_42_22","BA44_45","BA7","BA9"),
  CohenD = cohenF$fALFF
)

cortical_pos <- c("left lateral", "left medial")

pdf("revision/CohenD_fALFF1_GGSeg.pdf",height=5,width=5)
someData %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(cortical_pos),
             aes(fill = CohenD)) +
  scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  theme_void()
dev.off()


someData <- tibble(
  region = c("lateral occipital", "inferior temporal", "caudal anterior cingulate","postcentral","middle temporal","supramarginal",
  			"caudal middle frontal","superior temporal", "pars opercularis","superior parietal","rostral middle frontal"),
  ba =  c("BA17", "BA20_37","BA24","BA3_1_2_5","BA38","BA39_40","BA4_6","BA41_42_22","BA44_45","BA7","BA9"),
  CohenD = cohenR$ReHo
)

cortical_pos <- c("left lateral", "left medial")

pdf("revision/CohenD_ReHo_GGSeg.pdf",height=5,width=5)
someData %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(cortical_pos),
             aes(fill = CohenD)) +
  scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  theme_void()
dev.off()

