rm(list=ls())
suppressPackageStartupMessages({
library(here)
library(factoextra)
library(ade4)
library(tidyverse)
library(data.table)
library(FactoMineR)
library(ggpubr)
library(ggseg)
})

load(here("rawdata","Imaging_Normalized.RData"))

tmp <- Imaging %>%
		select(Region, fALFF1, ReHo1) %>%
      pivot_longer(!Region, names_to = "fMRI", values_to = "Value") %>%
      as.data.frame() %>%
      group_by(Region,fMRI) %>%
    dplyr::summarize(Mean = mean(Value, na.rm=TRUE)) %>%
    arrange(fMRI) %>%
    as.data.frame() 

someData <- tibble(
  region = rep(c("lateral occipital", "inferior temporal", "caudal anterior cingulate","postcentral","middle temporal","supramarginal",
  			"caudal middle frontal","superior temporal", "pars opercularis","superior parietal","rostral middle frontal"),2),
  ba =  rep(c("BA17", "BA20_37","BA24","BA3_1_2_5","BA38","BA39_40","BA4_6","BA41_42_22","BA44_45","BA7","BA9"),2),
  p = tmp$Mean,
  fMRI = tmp$fMRI
)

cortical_pos <- c("left lateral", "left medial")

pdf("revision/Distribution_Imaging_fALFF1_Seg.pdf",height=5,width=5)
someData %>%
  filter(fMRI == "fALFF1") %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(cortical_pos),
             aes(fill = p)) +
  scale_fill_gradientn(colours=c("#FFFFFFFF","#FF0000FF")) +
  theme_void()
dev.off()


pdf("revision/Distribution_Imaging_ReHo1_Seg.pdf",height=5,width=5)
someData %>%
  filter(fMRI == "ReHo1") %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(cortical_pos),
             aes(fill = p)) +
  scale_fill_gradientn(colours=c("#FFFFFFFF","#FF0000FF")) +
  theme_void()
dev.off()