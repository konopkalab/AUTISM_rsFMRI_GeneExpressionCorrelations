# Libraries and codes
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(made4))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(psych))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(matrixStats))
source("UTILS/Utils.R")

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

# Loading Imaging data
load(here("rawdata","Imaging_Normalized.RData"))

# Reorder the expression data
p_ctl1 <- p_ctl[,match(demoCTL$ID,colnames(p_ctl))]
p_asd1 <- p_asd[,match(demoASD$ID,colnames(p_asd))]

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

# Downsample analysis
resCTL_PanCort<- vector("list",length = B) 
resASD_PanCort<- vector("list",length = B)

#### Correlation analysis
for (i in 1:B)
  {
		resCTL_PanCort[[i]] <- FastCor(p_ctl1,as.numeric(tmpCTL_PanCort[[i]]$ReHo1),method="spearman",alternative="two.sided",cores=20,override=TRUE) %>%
      								      as.data.frame()	%>%
									          rownames_to_column('Gene')
	
		resASD_PanCort[[i]] <- FastCor(p_asd1,as.numeric(tmpASD_PanCort[[i]]$ReHo1),method="spearman",alternative="two.sided",cores=20,override=TRUE) %>%
      								      as.data.frame()	%>%
									          rownames_to_column('Gene') 
	}

#### Combine the data
CTL_ReHo<- resCTL_PanCort %>% 
            bind_cols() %>%
            #filter_at(vars(starts_with("Pval")), all_vars(. < 0.05)) %>%
            #filter_at(vars(starts_with("Rho")), all_vars(. < 0 | . > 0)) %>%
            mutate(Rho_CTL_ReHo = rowMeans(select(., starts_with("Rho")), na.rm = TRUE), Pval_CTL_ReHo = rowMeans(select(., starts_with("Pval")), na.rm = TRUE)) %>%
            mutate(FDR_CTL_ReHo = p.adjust(Pval_CTL_ReHo,method="BH")) %>%
            select("Gene","Rho_CTL_ReHo","Pval_CTL_ReHo","FDR_CTL_ReHo")


ASD_ReHo <- resASD_PanCort %>% 
            bind_cols() %>%
            mutate(Rho_ASD_ReHo = rowMeans(select(., starts_with("Rho")), na.rm = TRUE), Pval_ASD_ReHo = rowMeans(select(., starts_with("Pval")), na.rm = TRUE)) %>%
            mutate(FDR_ASD_ReHo = p.adjust(Pval_ASD_ReHo,method="BH")) %>%
            select("Gene","Rho_ASD_ReHo","Pval_ASD_ReHo","FDR_ASD_ReHo")

# Differential Correlation analysis             
DiffCor_ReHo <- left_join(CTL_ReHo,ASD_ReHo,by="Gene")%>%
                  mutate(DiffCor_ReHo_P = paired.r(Rho_CTL_ReHo,Rho_ASD_ReHo,NULL, 302, 360,twotailed=TRUE)$p,
                         DiffCor_ReHo_Z = paired.r(Rho_CTL_ReHo,Rho_ASD_ReHo,NULL, 302, 360,twotailed=TRUE)$z) 

dir.create("Permuted_Matches_Outputs/")
save(CTL_ReHo,ASD_ReHo,DiffCor_ReHo, file="Permuted_Matches_Outputs/Permuted_Matches_ReHo.RData")


#### Combine the data
CTL_ReHo<- resCTL_PanCort %>% 
            bind_cols() %>%
            filter_at(vars(starts_with("Pval")), all_vars(. < 0.05)) %>%
            filter_at(vars(starts_with("Rho")), all_vars(. < 0 | . > 0)) %>%
            mutate(Rho_CTL_ReHo = rowMeans(select(., starts_with("Rho")), na.rm = TRUE), Pval_CTL_ReHo = rowMeans(select(., starts_with("Pval")), na.rm = TRUE)) %>%
            mutate(FDR_CTL_ReHo = p.adjust(Pval_CTL_ReHo,method="BH")) %>%
            select("Gene","Rho_CTL_ReHo","Pval_CTL_ReHo","FDR_CTL_ReHo")


ASD_ReHo <- resASD_PanCort %>% 
            bind_cols() %>%
            mutate(Rho_ASD_ReHo = rowMeans(select(., starts_with("Rho")), na.rm = TRUE), Pval_ASD_ReHo = rowMeans(select(., starts_with("Pval")), na.rm = TRUE)) %>%
            mutate(FDR_ASD_ReHo = p.adjust(Pval_ASD_ReHo,method="BH")) %>%
            select("Gene","Rho_ASD_ReHo","Pval_ASD_ReHo","FDR_ASD_ReHo")

# Differential Correlation analysis             
DiffCor_ReHo_Specific <- left_join(CTL_ReHo,ASD_ReHo,by="Gene")%>%
                          mutate(DiffCor_ReHo_P = paired.r(Rho_CTL_ReHo,Rho_ASD_ReHo,NULL, 302, 360,twotailed=TRUE)$p,
                                 DiffCor_ReHo_Z = paired.r(Rho_CTL_ReHo,Rho_ASD_ReHo,NULL, 302, 360,twotailed=TRUE)$z) 

dir.create("Permuted_Matches_Outputs/")
save(DiffCor_ReHo_Specific, file="Permuted_Matches_Outputs/Permuted_Matches_ReHo_Specific.RData")

sessionInfo()






