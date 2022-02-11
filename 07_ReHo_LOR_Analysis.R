# Libraries and codes
suppressPackageStartupMessages({
library(reshape2)
library(purrr)
library(knitr)
library(broom)
library(tibble)
library(tidyverse)
library(here)
library(matrixStats)
library(preprocessCore)
library(openxlsx)
library(psych)
})
source("UTILS/Utils.R")

# Create Directory
dir.create("Permuted_Matches_LoRo_Outputs/")

# loading genomic data and formatting expression and demographic data
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

# Loading ReHo data
load(here("rawdata","Imaging_Normalized.RData"))

####################
### LOR Analysis ###
####################

# Create LOR demographic table
demoCTL$Region <- as.factor(demoCTL$Region)
N = unique(Imaging$Region) ## N region
tmpDemo_CTL <- vector("list", length = 11)
tmpDemo_ASD <- vector("list", length = 11)
tmpCTL_exp <- vector("list", length = 11)
tmpASD_exp <- vector("list", length = 11)
names_CTL <- vector("list", length = 11)
names_ASD <- vector("list", length = 11)

for (j in 1:11) 
      {
      tmpDemo_CTL[[j]] <- demoCTL %>%
                            filter(Region != paste(N[[j]]))%>%
                            droplevels()
      tmpDemo_ASD[[j]] <- demoASD %>%
                            filter(Region != paste(N[[j]]))%>%
                            droplevels()

      tmpCTL_exp[[j]] <- p_ctl[,colnames(p_ctl) %in% tmpDemo_CTL[[j]]$ID]

      tmpASD_exp[[j]] <- p_asd[,colnames(p_asd) %in% tmpDemo_ASD[[j]]$ID]

      names_CTL[[j]] <- sort(unique(tmpDemo_CTL[[j]]$Region)) 

      names_ASD[[j]] <- sort(unique(tmpDemo_ASD[[j]]$Region)) 

      }


# Nested loop for LOR matched ReHo values
B=200  ## select number of times for leave-multiple-out method * region
tmpCTL_PanCort <- vector("list", length = B)
for (j in 1:11)
  { 
    for (i in 1:B) 
      {

      tmpCTL_PanCort[[j]][[i]] <- Imaging %>%
                                  filter(Region %in% names_CTL[[j]]) %>%
                                  group_by(Region) %>%
                                  nest() %>%           
                                  arrange(Region, names_CTL[[j]]) %>% 
                                  mutate(n = sapply(split(tmpDemo_CTL[[j]],tmpDemo_CTL[[j]]$Region), nrow)) %>% 
                                  mutate(samp = map2(data, n, sample_n)) %>% 
                                  select(Region, samp) %>%
                                  unnest() %>%
                                  as.data.frame()
  }
}


# Nested loop for LOR matched ReHo values
B=200  ## select number of times for leave-multiple-out method * region
tmpASD_PanCort <- vector("list", length = B)
for (j in 1:11)
  { 
    for (i in 1:B) 
      {

      tmpASD_PanCort[[j]][[i]] <- Imaging %>%
                                  filter(Region %in% names_ASD[[j]]) %>%
                                  group_by(Region) %>%
                                  nest() %>%           
                                  arrange(Region, names_ASD[[j]]) %>% 
                                  mutate(n = sapply(split(tmpDemo_ASD[[j]],tmpDemo_ASD[[j]]$Region), nrow)) %>% 
                                  mutate(samp = map2(data, n, sample_n)) %>% 
                                  select(Region, samp) %>%
                                  unnest() %>%
                                  as.data.frame()
  }
}


#### Correlation analysis
resCTL_PanCort<- vector("list",length = 11) 
resASD_PanCort<- vector("list",length = 11)
for (i in 1:B)
{
  for (j in 1:11)
  {
    resCTL_PanCort[[j]][[i]] <- FastCor(tmpCTL_exp[[j]],as.numeric(tmpCTL_PanCort[[j]][[i]]$ReHo1),method="spearman",alternative="two.sided",cores=20,override=TRUE) %>%
                                  as.data.frame() %>%
                                  rownames_to_column('Gene')
  
    resASD_PanCort[[j]][[i]] <- FastCor(tmpASD_exp[[j]],as.numeric(tmpASD_PanCort[[j]][[i]]$ReHo1),method="spearman",alternative="two.sided",cores=20,override=TRUE) %>%
                                  as.data.frame() %>%
                                  rownames_to_column('Gene') 
  }
}
save(resCTL_PanCort,resASD_PanCort,file="Permuted_Matches_LoRo_Outputs/TemporaryFiles_LoRo_ReHo.RData")

# Combine data
CTL_ReHo_LoRo<- vector("list",length = 11) 
ASD_ReHo_LoRo<- vector("list",length = 11) 
DiffCor_ReHo_LoRo<- vector("list",length = 11) 

for (j in 1:11)
{
  CTL_ReHo_LoRo[[j]] <-  resCTL_PanCort[[j]] %>% 
                          bind_cols() %>%
                          mutate(Rho_CTL_ReHo = rowMeans(select(., starts_with("Rho")), na.rm = TRUE), Pval_CTL_ReHo = rowMeans(select(., starts_with("Pval")), na.rm = TRUE)) %>%
                          mutate(FDR_CTL_ReHo = p.adjust(Pval_CTL_ReHo,method="BH")) %>%
                          select("Gene","Rho_CTL_ReHo","Pval_CTL_ReHo","FDR_CTL_ReHo")
                          
  ASD_ReHo_LoRo[[j]] <- resASD_PanCort[[j]] %>% 
                          bind_cols() %>%
                          mutate(Rho_ASD_ReHo = rowMeans(select(., starts_with("Rho")), na.rm = TRUE), Pval_ASD_ReHo = rowMeans(select(., starts_with("Pval")), na.rm = TRUE)) %>%
                          mutate(FDR_ASD_ReHo = p.adjust(Pval_ASD_ReHo,method="BH")) %>%
                          select("Gene","Rho_ASD_ReHo","Pval_ASD_ReHo","FDR_ASD_ReHo")

  DiffCor_ReHo_LoRo[[j]] <- left_join(CTL_ReHo_LoRo[[j]],ASD_ReHo_LoRo[[j]],by="Gene")%>%
                            mutate(DiffCor_ReHo_P = paired.r(Rho_CTL_ReHo,Rho_ASD_ReHo,NULL, 302, 360,twotailed=TRUE)$p,
                                   DiffCor_ReHo_Z = paired.r(Rho_CTL_ReHo,Rho_ASD_ReHo,NULL, 302, 360,twotailed=TRUE)$z)    
}

names(CTL_ReHo_LoRo) <- paste("LoRo_",N,sep="")
names(ASD_ReHo_LoRo) <- paste("LoRo_",N,sep="")
names(DiffCor_ReHo_LoRo) <- paste("LoRo_",N,sep="")

# save lists
save(CTL_ReHo_LoRo,ASD_ReHo_LoRo,DiffCor_ReHo_LoRo, file="Permuted_Matches_LoRo_Outputs/Permuted_Matches_LoRo_ReHo.RData")


