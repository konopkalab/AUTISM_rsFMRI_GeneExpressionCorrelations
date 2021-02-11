rm(list=ls())
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggjoy))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(psych))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(VennDiagram))

# Create Output directories
folder_names <- c("integrative_results")
sapply(folder_names, dir.create)

# load data you need for intersection
load(here("Permuted_Matches_Outputs", "Permuted_Matches_fALFF.RData"))
load(here("Permuted_Matches_Outputs", "Permuted_Matches_ReHo.RData"))

# Raneming columns
DiffCor_fALFF <- DiffCor_fALFF %>%
                    mutate(Rsq_CTL_fALFF = Rho_CTL_fALFF^2,
                           Rsq_ASD_fALFF = Rho_ASD_fALFF^2)

DiffCor_ReHo <- DiffCor_ReHo %>%
                  mutate(Rsq_CTL_ReHo = Rho_CTL_ReHo^2,
                         Rsq_ASD_ReHo = Rho_ASD_ReHo^2)

# Combine all data
DiffCor_Combined <- Reduce(dplyr::full_join, list(DiffCor_fALFF, DiffCor_ReHo))
openxlsx::write.xlsx(DiffCor_Combined, file = "integrative_results/DiffCor_Combined_TableS2.xlsx", colNames = TRUE, borders = "columns",sheetName="Full Table")
save(DiffCor_fALFF,DiffCor_ReHo,DiffCor_Combined,file = "integrative_results/DiffCor_Tables.RData")

# CTL specific only
Only_CTL_fALFF <- DiffCor_fALFF %>% 
                        filter(FDR_CTL_fALFF < 0.05) %>%
                        as.data.frame()

Only_CTL_ReHo <- DiffCor_ReHo %>% 
                        filter(FDR_CTL_ReHo < 0.05) %>%
                        as.data.frame()

save(Only_CTL_fALFF,Only_CTL_ReHo,file = "integrative_results/Only_CTL_Significant.RData")


# ASD specific only
Only_ASD_fALFF <- DiffCor_fALFF %>% 
                        filter(FDR_ASD_fALFF < 0.05) %>%
                        as.data.frame()

Only_ASD_ReHo <- DiffCor_ReHo %>% 
                        filter(FDR_ASD_ReHo < 0.05) %>%
                        as.data.frame()

save(Only_ASD_fALFF,Only_ASD_ReHo,file = "integrative_results/Only_ASD_Significant.RData")


# Diff Cor specific
DiffCor_Sign_fALFF <- DiffCor_fALFF %>% 
                        filter(FDR_CTL_fALFF < 0.05 & DiffCor_fALFF_P < 0.05) %>%
                        as.data.frame()

DiffCor_Sign_ReHo <- DiffCor_ReHo %>% 
                        filter(FDR_CTL_ReHo < 0.05 & DiffCor_ReHo_P < 0.05) %>%
                        as.data.frame()

save(DiffCor_Sign_fALFF,DiffCor_Sign_ReHo,file = "integrative_results/DiffCor_Significant.RData")

xlsx::write.xlsx(DiffCor_Sign_fALFF, file="integrative_results/DiffCor_Significant_TableS3.xlsx",sheetName = "fALFF",row.names=FALSE, showNA=FALSE)
xlsx::write.xlsx(DiffCor_Sign_ReHo, file="integrative_results/DiffCor_Significant_TableS3.xlsx",sheetName = "ReHo",row.names=FALSE, showNA=FALSE,append=TRUE)

