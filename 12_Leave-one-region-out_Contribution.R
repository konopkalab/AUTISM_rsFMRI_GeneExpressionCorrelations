rm(list=ls())
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(ade4))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(FactoMineR))
suppressPackageStartupMessages(library(ggpubr))

load(here("Permuted_Matches_LoRo_Outputs", "Permuted_Matches_LoRo_fALFF.RData"))
genes <- read.table(here("integrative_results","ASD_fMRI_Genes.txt"),header=T,sep="\t")

genes <- genes %>%
          filter(Direction == "AllGene")


new_names <- gsub("LoRo_","",names(DiffCor_fALFF_LoRo))

# Define colors
paired <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00","#CAB2D6","#6A3D9A","#FF99F8")

# Visualizations
pdf("integrative_visualizations/PCA_LoR_Contribution_fALFF.pdf", width=5, height=4)
map(DiffCor_fALFF_LoRo, ~ (.x %>% select(Gene,DiffCor_fALFF_Z))) %>% 
      map2(new_names, ~setnames(.x, 'DiffCor_fALFF_Z', .y)) %>%
      reduce(left_join, by = "Gene") %>% 
      filter(Gene %in% genes$Gene) %>% 
      column_to_rownames("Gene") %>%
      t() %>%
      prcomp(scale = FALSE) %>%
      fviz_pca_ind(
             col.ind = "contrib",
             repel = TRUE,
             alpha.ind="contrib"
             ) +
      scale_color_gradient2(low="grey", mid="steelblue",high="red", midpoint=15)+
             xlim(-15,15) +
             ylim(-10,10)
dev.off()

# Save into an excel file
map(DiffCor_fALFF_LoRo, ~ (.x %>% select(Gene,DiffCor_fALFF_Z))) %>% 
      map2(new_names, ~setnames(.x, 'DiffCor_fALFF_Z', .y)) %>%
      reduce(left_join, by = "Gene") %>% 
      filter(Gene %in% genes$Gene) %>% 
      column_to_rownames("Gene") %>%
      openxlsx::write.xlsx(file = "supp_tables/DiffCor_fALFF_LoRo_Values.xlsx", colNames = TRUE, borders = "columns")

pdf("integrative_visualizations/Boxplot_LoR_Contribution_fALFF.pdf", width=5, height=3)
 map(DiffCor_fALFF_LoRo, ~ (.x %>% select(Gene,DiffCor_fALFF_Z))) %>% 
      map2(new_names, ~setnames(.x, 'DiffCor_fALFF_Z', .y)) %>%
       reduce(left_join, by = "Gene") %>% 
       filter(Gene %in% genes$Gene) %>% 
       column_to_rownames("Gene") %>%
       melt() %>%
       rename(Region = variable, Zscore = value) %>%
       mutate(absZ = abs(Zscore)) %>%
       ggviolin(
         x = "Region", 
         y = "Zscore", 
         color = "Region", 
           add = c("mean_sd"), 
           legend = "none",
           palette = paired) +
     rotate_x_text(angle = 45)+
          stat_compare_means(label = "p.signif", 
                          method = "wilcox.test",
                          method.args = list(alternative = "less"),
                          ref.group = ".all.")+
    xlab("") 
dev.off()

pdf("integrative_visualizations/PCA_LoR_Contribution_fALFF_BarPlot.pdf", width=5, height=4)
map(DiffCor_fALFF_LoRo, ~ (.x %>% select(Gene,DiffCor_fALFF_Z))) %>% 
      map2(new_names, ~setnames(.x, 'DiffCor_fALFF_Z', .y)) %>%
      reduce(left_join, by = "Gene") %>% 
      filter(Gene %in% genes$Gene) %>% 
      column_to_rownames("Gene") %>%
      t() %>%
      prcomp(scale = FALSE) %>%
      fviz_contrib(choice = "ind",axes = 1:2, top = 12,fill = "white", color = "black")+
      theme_classic()+
      rotate_x_text(angle = 45)
dev.off()


# Gain and lost
lost <- map(DiffCor_fALFF_LoRo, ~ (.x %>% filter(Pval_CTL_fALFF < 0.05, DiffCor_fALFF_P < 0.05, !Gene %in% genes$Gene))) 
sign <- map(DiffCor_fALFF_LoRo, ~ (.x %>% filter(Pval_CTL_fALFF < 0.05, DiffCor_fALFF_P < 0.05)))

fALFF_Gained <- vector("list", length = 11)
for (j in 1:11) 
      {
      fALFF_Gained[[j]] <- genes %>%
                            filter(!(Gene %in% sign[[j]]$Gene))
      }

save(fALFF_Gained, file="Permuted_Matches_LoRo_Outputs/fALFF_Gained.RData")

pdf("integrative_visualizations/PCA_LoR_Gained_fALFF.pdf", width=4, height=3)
data.frame(Region = new_names, Gain = sapply(fALFF_Gained,nrow), Lost = sapply(lost,nrow)) %>%
            #mutate(Scale = as.numeric(scale(Gain))) %>%
            arrange(Region) %>%
            select(Region,Gain) %>%
            melt() %>%
            ggbarplot(x = "Region", 
                              y = "value",
                        color = "black",            
                        #palette = "jco",            
                        #sort.val = "asc",          
                        sort.by.groups = FALSE,     
                        x.text.angle = 45,          
                        ylab = "# of Genes Gained",
                        rotate = FALSE,
                          ggtheme = theme_minimal()
                      )+
            #ylim(3,-3)+
            theme(legend.position="none")+
            xlab("")+
            ylim(0,400)
dev.off()


# Loading for ReHo
rm(list=ls())
load(here("Permuted_Matches_LoRo_Outputs", "Permuted_Matches_LoRo_ReHo.RData"))
genes <- read.table(here("integrative_results","ASD_fMRI_Genes.txt"),header=T,sep="\t")

genes <- genes %>%
          filter(Direction == "AllGene")


new_names <- gsub("LoRo_","",names(DiffCor_ReHo_LoRo))

pdf("integrative_visualizations/PCA_LoR_Contribution_ReHo.pdf", width=5, height=4)
map(DiffCor_ReHo_LoRo, ~ (.x %>% select(Gene,DiffCor_ReHo_Z))) %>% 
      map2(new_names, ~setnames(.x, 'DiffCor_ReHo_Z', .y)) %>%
      reduce(left_join, by = "Gene") %>% 
      filter(Gene %in% genes$Gene) %>% 
      column_to_rownames("Gene") %>%
      t() %>%
      prcomp(scale = FALSE) %>%
      fviz_pca_ind(
             col.ind = "contrib",
             repel = TRUE,
             alpha.ind="contrib"
             ) +
      scale_color_gradient2(low="grey", mid="steelblue",high="red", midpoint=8)+
             xlim(-15,15) +
             ylim(-10,10)
dev.off()

# Save into an excel file
map(DiffCor_ReHo_LoRo, ~ (.x %>% select(Gene,DiffCor_ReHo_Z))) %>% 
      map2(new_names, ~setnames(.x, 'DiffCor_ReHo_Z', .y)) %>%
      reduce(left_join, by = "Gene") %>% 
      filter(Gene %in% genes$Gene) %>% 
      column_to_rownames("Gene") %>%
      openxlsx::write.xlsx(file = "supp_tables/DiffCor_ReHo_LoRo_Values.xlsx", colNames = TRUE, borders = "columns")

paired <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00","#CAB2D6","#6A3D9A","#FF99F8")

pdf("integrative_visualizations/Boxplot_LoR_Contribution_ReHo.pdf", width=5, height=3)
 map(DiffCor_ReHo_LoRo, ~ (.x %>% select(Gene,DiffCor_ReHo_Z))) %>% 
      map2(new_names, ~setnames(.x, 'DiffCor_ReHo_Z', .y)) %>%
      reduce(left_join, by = "Gene") %>% 
      filter(Gene %in% genes$Gene) %>% 
      column_to_rownames("Gene")  %>%
       melt() %>%
       rename(Region = variable, Zscore = value) %>%
       mutate(absZ = abs(Zscore)) %>%
       ggviolin(
         x = "Region", 
         y = "Zscore", 
         color = "Region", 
           add = c("mean_sd"), 
           legend = "none",
           palette = paired) +
     rotate_x_text(angle = 45)+
          stat_compare_means(label = "p.signif", 
                          method = "wilcox.test",
                          method.args = list(alternative = "less"),
                          ref.group = ".all.")+
    xlab("") 
dev.off()

pdf("integrative_visualizations/PCA_LoR_Contribution_ReHo_BarPlot.pdf", width=5, height=4)
map(DiffCor_ReHo_LoRo, ~ (.x %>% select(Gene,DiffCor_ReHo_Z))) %>% 
      map2(new_names, ~setnames(.x, 'DiffCor_ReHo_Z', .y)) %>%
      reduce(left_join, by = "Gene") %>% 
      filter(Gene %in% genes$Gene) %>% 
      column_to_rownames("Gene") %>%
      t() %>%
      prcomp(scale = FALSE) %>%
      fviz_contrib(choice = "ind",axes = 1:2, top = 12,fill = "white", color = "black")+
      theme_classic()+
      rotate_x_text(angle = 45)
dev.off()


# Gain and lost
lost <- map(DiffCor_ReHo_LoRo, ~ (.x %>% filter(Pval_CTL_ReHo < 0.05, DiffCor_ReHo_P < 0.05, !Gene %in% genes$Gene))) 
sign <- map(DiffCor_ReHo_LoRo, ~ (.x %>% filter(Pval_CTL_ReHo < 0.05, DiffCor_ReHo_P < 0.05)))

ReHo_Gained <- vector("list", length = 11)
for (j in 1:11) 
      {
      ReHo_Gained[[j]] <- genes %>%
                            filter(!(Gene %in% sign[[j]]$Gene))
      }
save(ReHo_Gained, file="Permuted_Matches_LoRo_Outputs/ReHo_Gained.RData")


pdf("integrative_visualizations/PCA_LoR_Gained_ReHo.pdf", width=4, height=3)
data.frame(Region = new_names, Gain = sapply(ReHo_Gained,nrow), Lost = sapply(lost,nrow)) %>%
            #mutate(Scale = as.numeric(scale(Gain))) %>%
            arrange(Region) %>%
            select(Region,Gain) %>%
            melt() %>%
            ggbarplot(x = "Region", 
                              y = "value",
                        color = "black",            
                        #palette = "jco",            
                        #sort.val = "asc",          
                        sort.by.groups = FALSE,     
                        x.text.angle = 45,          
                        ylab = "# of Genes Gained",
                        rotate = FALSE,
                          ggtheme = theme_minimal()
                      )+
            #ylim(3,-3)+
            theme(legend.position="none")+
            xlab("")+
            ylim(0,400)
dev.off()
