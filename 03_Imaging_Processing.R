# Covariate plot
rm(list=ls())
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(matrixStats))

###################################### 
# Loading fALFF data and normalizing #
######################################
fALFF1 <- read.table(here("rawdata","ABIDEIandIIChosenScansBrodmannValues_equalized_global_fALFF.csv"),sep=",",header=T)
colnames(fALFF1) <- gsub("MeanRegion","BA",colnames(fALFF1))

# Remove un-used columns
fALFF1 <- fALFF1 %>% 
      mutate(Diagnosis = case_when(DX == 1 ~ "ASD", DX == 2 ~ "CTL"),
             Sex = case_when(SEX == 2 ~ "F", SEX == 1 ~ "M")) %>%
      select(-DX, -SEX) %>%
      select(Sex, everything()) %>%
      select(Diagnosis, everything()) %>%
      select(ID, everything()) %>%
#      filter(AGE_AT_SCAN > 18) %>%
      as.data.frame() 

openxlsx::write.xlsx(fALFF1, file = "supp_tables/fALFF1_ABIDE_Values.xlsx", colNames = TRUE, borders = "columns")

# Mutate columns as factor
cols <- c("ID","Diagnosis","Sex","ABIDE","SITE")
fALFF1 %<>% mutate_at(cols, funs(factor(.)))


vox <- read.table(here("rawdata","Regions_NumberOfVoxels.txt"),header=T,sep="\t")

# Adjust fALFF data
tmpA <- as.data.frame(t(fALFF1[,7:17]/rowWeightedMeans(as.matrix(fALFF1[,7:17]),w=vox$Voxels)))
pd_sva <- droplevels(fALFF1[c(3,5,6)])
avebeta.lm<-lapply(1:nrow(tmpA), function(x)
          {
              lm(unlist(tmpA[x,])~., data = pd_sva)
          })
residuals<-lapply(avebeta.lm, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
tmpA_adj<-residuals+matrix(apply(tmpA, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(tmpA_adj)<-rownames(tmpA)
fALFF1 <- cbind(fALFF1[,1:6],t(tmpA_adj))

# Filtering, reshaping and matching with expression data
fALFF1_filt <- droplevels(fALFF1[c(1,2,7:17)])

fALFF1_Data <- fALFF1_filt %>% 
         gather(Region,fALFF1,-ID, -Diagnosis) %>%
         unite_("NewID", c("ID","Region"),remove = FALSE)

# Filter for same subject - same region 
fALFF_CTL <- fALFF1_Data %>% 
              filter(Diagnosis == "CTL") %>%
              droplevels() %>%
              arrange(Region)

save(fALFF1_Data, fALFF_CTL, file = "rawdata/fALFF_CTL_Shaped.RData")

##################################### 
# Loading ReHo data and normalizing #
#####################################
ReHo1 <- read.table(here("rawdata","ABIDEIandIIChosenScansBrodmannValues_equalized_global_ReHo.csv"),sep=",",header=T)
colnames(ReHo1) <- gsub("MeanRegion","BA",colnames(ReHo1))

# Remove un-used columns
ReHo1 <- ReHo1 %>% 
      mutate(Diagnosis = case_when(DX == 1 ~ "ASD", DX == 2 ~ "CTL"),
             Sex = case_when(SEX == 2 ~ "F", SEX == 1 ~ "M")) %>%
      select(-DX, -SEX) %>%
      select(Sex, everything()) %>%
      select(Diagnosis, everything()) %>%
      select(ID, everything()) %>%
#      filter(AGE_AT_SCAN > 18) %>%
      as.data.frame() 
openxlsx::write.xlsx(ReHo1, file = "supp_tables/ReHo1_ABIDE_Values.xlsx", colNames = TRUE, borders = "columns")

# Mutate columns as factor
cols <- c("ID","Diagnosis","Sex","ABIDE","SITE")
ReHo1 %<>% mutate_at(cols, funs(factor(.)))


vox <- read.table(here("rawdata","Regions_NumberOfVoxels.txt"),header=T,sep="\t")

# Adjust fALFF data
tmpA <- as.data.frame(t(ReHo1[,7:17]/rowWeightedMeans(as.matrix(ReHo1[,7:17]),w=vox$Voxels)))
pd_sva <- droplevels(ReHo1[c(3,5,6)])
avebeta.lm<-lapply(1:nrow(tmpA), function(x)
          {
              lm(unlist(tmpA[x,])~., data = pd_sva)
          })
residuals<-lapply(avebeta.lm, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
tmpA_adj<-residuals+matrix(apply(tmpA, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(tmpA_adj)<-rownames(tmpA)
ReHo1 <- cbind(ReHo1[,1:6],t(tmpA_adj))

# Filtering, reshaping and matching with expression data
ReHo1_filt <- droplevels(ReHo1[c(1,2,7:17)])

ReHo1_Data <- ReHo1_filt %>% 
         gather(Region,ReHo1,-ID, -Diagnosis) %>%
         unite_("NewID", c("ID","Region"),remove = FALSE)

# Filter for same subject - same region 
ReHo_CTL <- ReHo1_Data %>% 
              filter(Diagnosis == "CTL") %>%
              droplevels()%>%
              arrange(Region)
              
save(ReHo1_Data,ReHo_CTL, file = "rawdata/ReHo_CTL_Shaped.RData")

########################
### Joing the tables ###
########################
Imaging <- Reduce(dplyr::full_join, list(fALFF_CTL,ReHo_CTL))%>%
           na.omit() # removing subjects specific
save(Imaging, file = "rawdata/Imaging_Normalized.RData")

Imaging_ASD <- Reduce(dplyr::full_join, list(fALFF1_Data,ReHo1_Data))%>%
           na.omit() # removing subjects specific
save(Imaging_ASD, file = "rawdata/Imaging_with_ASD_Normalized.RData")

sessionInfo()