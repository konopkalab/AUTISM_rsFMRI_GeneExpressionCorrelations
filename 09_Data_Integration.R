rm(list=ls())
suppressPackageStartupMessages({
library(tidyverse)
library(here)
library(scran)
})

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

# Combine all data and combine P-value
DiffCor_Combined <- Reduce(dplyr::full_join, list(DiffCor_fALFF, DiffCor_ReHo)) %>% 
                    mutate(Pval_CTL_Comb = combinePValues(Pval_CTL_fALFF,Pval_CTL_ReHo),
                    	   Pval_ASD_Comb = combinePValues(Pval_ASD_fALFF,Pval_ASD_ReHo),
                    	   DiffCor_Comb_P = combinePValues(DiffCor_fALFF_P,DiffCor_ReHo_P)) %>%
                    mutate(FDR_CTL_Comb = p.adjust(Pval_CTL_Comb,method="BH"),
                    	   FDR_ASD_Comb = p.adjust(Pval_ASD_Comb,method="BH")) %>%
                    as.data.frame()

openxlsx::write.xlsx(DiffCor_Combined, file = "integrative_results/DiffCor_Combined_TableS2.xlsx", 
                    colNames = TRUE, borders = "columns",sheetName="Full Table",overwrite=TRUE)

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
DiffCor_Sign_Final <- DiffCor_Combined %>% 
                        filter(FDR_CTL_Comb < 0.05 & DiffCor_Comb_P < 0.01 & sign(Rho_CTL_fALFF) == sign(Rho_CTL_ReHo)) %>% 
                        select(-Pval_CTL_fALFF,-FDR_CTL_fALFF,-Pval_ASD_fALFF,-FDR_ASD_fALFF,
                        	   -Pval_CTL_ReHo,-FDR_CTL_ReHo,-Pval_ASD_ReHo,-FDR_ASD_ReHo,
                        	   -DiffCor_fALFF_P,-DiffCor_ReHo_P) %>%
                        mutate(Direction= case_when(Rho_CTL_fALFF > 0 ~ "Positive", Rho_CTL_fALFF < 0  ~ "Negative")) %>%
                        as.data.frame()

CTL_Sign_Final <- DiffCor_Combined %>% 
                        filter(FDR_CTL_Comb < 0.05 & sign(Rho_CTL_fALFF) == sign(Rho_CTL_ReHo)) %>% 
                        select(-Pval_CTL_fALFF,-FDR_CTL_fALFF,-Pval_ASD_fALFF,-FDR_ASD_fALFF,
                               -Pval_CTL_ReHo,-FDR_CTL_ReHo,-Pval_ASD_ReHo,-FDR_ASD_ReHo,
                               -DiffCor_fALFF_P,-DiffCor_ReHo_P) %>%
                        mutate(Direction= case_when(Rho_CTL_fALFF > 0 ~ "Positive", Rho_CTL_fALFF < 0  ~ "Negative")) %>%
                        as.data.frame()

save(CTL_Sign_Final,DiffCor_Sign_Final,file = "integrative_results/DiffCor_Significant.RData")

openxlsx::write.xlsx(DiffCor_Sign_Final, file = "integrative_results/DiffCor_Significant_TableS3.xlsx", 
                    colNames = TRUE, borders = "columns",sheetName="Full Table",overwrite=TRUE)

openxlsx::write.xlsx(CTL_Sign_Final, file = "integrative_results/CTL_Significant_TableS3.xlsx", 
                    colNames = TRUE, borders = "columns",sheetName="Full Table",overwrite=TRUE)
