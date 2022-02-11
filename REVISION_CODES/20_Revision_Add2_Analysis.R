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

load(here("Permuted_Matches_LoRo_Outputs", "Permuted_Matches_LoRo_fALFF.RData"))
load(here("Permuted_Matches_LoRo_Outputs", "Permuted_Matches_LoRo_ReHo.RData"))

genes <- read.table(here("integrative_results","ASD_fMRI_Genes.txt"),header=T,sep="\t")

genes <- genes %>%
          filter(Direction == "AllGene")


new_names <- gsub("LoRo_","",names(DiffCor_fALFF_LoRo))


fALFF <- map(DiffCor_fALFF_LoRo, ~ (.x %>% select(Gene,DiffCor_fALFF_Z))) %>% 
      map2(new_names, ~setnames(.x, 'DiffCor_fALFF_Z', .y)) %>%
      reduce(left_join, by = "Gene") %>% 
      filter(Gene %in% genes$Gene) %>% 
      pivot_longer(!Gene, names_to = "Region", values_to = "Zscore_fALFF") %>%
      as.data.frame()

ReHo <- map(DiffCor_ReHo_LoRo, ~ (.x %>% select(Gene,DiffCor_ReHo_Z))) %>% 
      map2(new_names, ~setnames(.x, 'DiffCor_ReHo_Z', .y)) %>%
      reduce(left_join, by = "Gene") %>% 
      filter(Gene %in% genes$Gene) %>% 
      pivot_longer(!Gene, names_to = "Region", values_to = "Zscore_ReHo") %>%
      as.data.frame()

LoRo_Combined <- Reduce(dplyr::full_join, list(fALFF, ReHo)) %>% 
				 rowwise() %>%
				 mutate(DiffCor_LoRo_Mean = mean(c(Zscore_fALFF,Zscore_ReHo))) %>%
				 select(Gene,Region,DiffCor_LoRo_Mean)


pdf("revision/Boxplot_LoRo_Contribution_Combined.pdf", width=5, height=3)
paired <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00","#CAB2D6","#6A3D9A","#FF99F8")
ggviolin(LoRo_Combined,
         x = "Region", 
         y = "DiffCor_LoRo_Mean", 
         color = "Region", 
           add = c("mean_sd"), 
           legend = "none",
           palette = paired) +
     rotate_x_text(angle = 45)+
          stat_compare_means(label = "p.format", 
                          method = "wilcox.test",
                          method.args = list(alternative = "less"),
                          ref.group = ".all.",
                          hide.ns = TRUE)
dev.off()


## GGseg
someData <- tibble(
  region = c("lateral occipital", "inferior temporal", "caudal anterior cingulate","postcentral","middle temporal","supramarginal",
  			"caudal middle frontal","superior temporal", "pars opercularis","superior parietal","rostral middle frontal"),
  ba =  c("BA17", "BA20_37","BA24","BA3_1_2_5","BA38","BA39_40","BA4_6","BA41_42_22","BA44_45","BA7","BA9"),
  p = c(0.00001,0.00001,1,1,0.0001,0.00001,0.00001,1,1,1,1)
)

pdf("revision/GGseg_LoRo.pdf",height=5,width=5)
someData %>%
  mutate(log = -log10(p)) %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),
             aes(fill = log)) +
  scale_fill_viridis_c(option = "cividis", direction = 1) +
  theme_void() 
dev.off()


pairwise <- radiant.basics::compare_means(
  LoRo_Combined %>% as.data.frame(),
  "Region",
  "DiffCor_LoRo_Mean",
  samples = "independent",
  alternative = "less",
  conf_lev = 0.95,
  comb = "",
  adjust = "none",
  test = "wilcox"
)$res %>%
mutate(log = -log10(p.value)) %>%
select(group1, group2, log) %>%
pivot_wider(names_from = group2, values_from = log) %>%
as.data.frame() %>%
column_to_rownames("group1")

pairwise[is.na(pairwise)] <- 1

cols = colorRampPalette(c("white", "red"))(30)

pdf("revision/LoRo_Contribution_Combined_PairwiseComparison.pdf",height=5,width=5)
pheatmap::pheatmap(pairwise, color = cols)
dev.off()
