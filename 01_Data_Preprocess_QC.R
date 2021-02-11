# Libraries and codes
suppressWarnings(suppressPackageStartupMessages(library(biomaRt)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(AnnotationDbi)))
suppressWarnings(suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86)))
suppressWarnings(suppressPackageStartupMessages(library(magrittr)))
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(here)))
suppressWarnings(suppressPackageStartupMessages(library(sva)))
suppressWarnings(suppressPackageStartupMessages(library(limma)))
suppressWarnings(suppressPackageStartupMessages(library(preprocessCore)))
suppressWarnings(suppressPackageStartupMessages(library(doParallel)))
suppressWarnings(suppressPackageStartupMessages(library(xlsx)))
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(purrr)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(gdata)))
suppressWarnings(suppressPackageStartupMessages(library(openxlsx)))
suppressWarnings(suppressPackageStartupMessages(library(DT)))
suppressWarnings(suppressPackageStartupMessages(library(future.apply)))


plan("multiprocess", workers = 3)


# Load expression data
load(here("rawdata","01_15_A-Dataset-Complete_08142018.RData"))

# Create a subject_region ID
datMeta$ID <- as.factor(paste(datMeta$subject, datMeta$region,sep="_"))
datSeq$ID <- as.factor(paste(datMeta$subject, datMeta$region,sep="_"))

# Select sample duplicated with highest RIN and remove remaining technical duplicates. 
datMeta <- datMeta %>% 
            dplyr::group_by(ID) %>% 
            dplyr::filter(RIN == max(RIN)) %>%
            as.data.frame()

datMeta <- datMeta[!(duplicated(datMeta$ID)),]
datSeq <- datSeq[rownames(datSeq)%in%datMeta$sample_id,]
datSeq <- datSeq[match(datMeta$sample_id,rownames(datSeq)),]

# Format Gene Symbol
exp <- rsem_gene
exp$Ensembl_ID <- gsub("\\..*","",rownames(exp))

# Annotate by gene names and select only protein coding
exp$Symbol <- mapIds(EnsDb.Hsapiens.v86,keys=exp$Ensembl_ID, column="GENENAME",keytype="GENEID",multiVals="first") # Get gene symbol
exp$BioType <- mapIds(EnsDb.Hsapiens.v86,keys=exp$Ensembl_ID, column="TXBIOTYPE",keytype="GENEID",multiVals="first") # Get biotype
exp <- exp[exp$BioType == "protein_coding",] # Select only protein coding genes
exp$BioType <- NULL
exp$Symbol[is.na(exp$Symbol)] <- exp$Ensembl_ID[is.na(exp$Symbol)] 
exp <- na.omit(exp)
agg <- aggregate(x = exp[, 1:854], by = list(ID = exp$Symbol), FUN = "max", na.rm = T) # Aggregate for multiple genes by mean
rownames(agg) <- agg$ID
agg$ID <- NULL

#
detach("package:ensembldb_2.8.0", unload=TRUE)

# Filter and change names
exp <- agg[,colnames(agg)%in%datMeta$sample_id,]
exp <- exp[,match(datMeta$sample_id,colnames(exp))]
colnames(exp) <- datMeta$ID
rownames(datMeta) <- datMeta$ID
rownames(datSeq) <- datSeq$ID

# Remove low expressed genes and log2 scale 
CPM <- future_apply(exp, 2, function(x) x/sum(as.numeric(x)) * 10^6)
exp <- as.data.frame(log2(CPM+1)) ## counts per million and log2 scaled data
pres = future_apply(exp>0.5,1,sum) 
idx = which(pres > 0.8*dim(exp)[2]) ## exp > 0.5 in 80% of samples
exp = exp[idx,]

# Using only ASD and CTL and filter for only region matching ABIDE data
pd <- datMeta %>%
        rownames_to_column("Temp") %>%
        filter(!(Diagnosis_ %in% c("Dup15q","NCTL"))) %>%
        filter(region %in% c("BA3-1-2-5","BA4-6","BA7","BA9","BA17","BA20-37","BA24","BA38","BA39-40","BA41-42-22","BA44-45")) %>%
        column_to_rownames("Temp") %>%
        droplevels()

pdSeq <- datSeq[rownames(datSeq)%in%rownames(pd),] #Filter datSeq
pdSeq$ID <- NULL
pd$ID <- NULL

# Filter the expression
dropInd <- which(names(exp) %in% rownames(pd))
exp <- exp[,dropInd]

# Covariate for regression
tmp <- as.data.frame(cbind(pd,pdSeq))
demo <- tmp %>% 
            dplyr::select(
                    Diagnosis_,
                    subject,
                    region,
                    seq_batch,
                    Brain_Bank_Source, 
                    Sex,  
                    Age, 
                    PMI, 
                    Ancestry_Genotype,
                    Previously_reported_RIN_CTX, 
                    picard_gcbias.AT_DROPOUT, 
                    star.deletion_length, 
                    picard_rnaseq.PCT_INTERGENIC_BASES, 
                    picard_insert.MEDIAN_INSERT_SIZE, 
                    picard_alignment.PCT_CHIMERAS, 
                    picard_alignment.PCT_PF_READS_ALIGNED, 
                    star.multimapped_percent, 
                    picard_rnaseq.MEDIAN_5PRIME_BIAS, 
                    star.unmapped_other_percent, 
                    picard_rnaseq.PCT_USABLE_BASES, 
                    star.uniquely_mapped_percent) %>%
            as.data.frame()

colnames(demo)[1:10] <- c("Diagnosis","Subject","Region","Batch","BrainBank","Sex","Age","PMI","Ancestry","RIN")
demo$PMI[is.na(demo$PMI)] <- median(demo$PMI, na.rm = T) # PMI imputation
demo$RIN[is.na(demo$RIN)] <- median(demo$RIN, na.rm = T) # RIN imputation
cols <- c("Diagnosis","Subject","Batch","BrainBank","Sex","Ancestry")
demo %<>% mutate_at(cols, funs(factor(.))) # Factorized for linear regression
rownames(demo) <- rownames(pd)
demo <- droplevels(demo)

# Adjust Expression with linear modeling
pd_sva <- demo %>%
            dplyr::select(-Diagnosis, -Region,-Subject) %>% #Removing ID, Region, Diagnosis from the model
            droplevels()

betas<-future_lapply(1:nrow(exp), function(x)
            {
                lm(unlist(exp[x,])~., data = pd_sva)
            })

residuals<-future_lapply(betas, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
expAdj <- residuals+matrix(future_apply(exp, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(expAdj)<-rownames(exp)

# add datSeq and datMeta
int <- intersect(datMeta$ID,rownames(demo))
int2 <- intersect(datSeq$ID,rownames(demo))

datMeta <- datMeta[datMeta$ID %in% int,]
datSeq <- datSeq[datSeq$ID %in% int2,]

datMeta<-datMeta[match(rownames(demo),datMeta$ID),]
datSeq<-datSeq[match(rownames(demo),datSeq$ID),]

save(exp,expAdj,demo,datMeta,datSeq,file="rawdata/Expression_Data.RData")

######################
#### DGE analysis ####
######################
folder_names <- c("dge")
sapply(folder_names, dir.create)

# load the data
load(here("rawdata","Expression_Data.RData"))
dat <- exp #Expression values unadjusted
pd <- demo #Demographic data sequencing covariates

# Check the structure of demographic 
str(demo)

# Running DGE
pd$Subject <- NULL
#pd$Region <- NULL
pd <- droplevels(pd)

total <- as.formula(paste("~", paste(names(pd),collapse="+")))
mod <- model.matrix(total, data =pd)
fitLM = lmFit(dat,mod,method="robust")
fitEb = eBayes(fitLM)
DGE = topTable(fitEb, coef = "DiagnosisASD",number=nrow(dat));
sign <- DGE[DGE$adj.P.Val < 0.05 & abs(DGE$logFC) > 0.3,]
xlsx::write.xlsx(DGE, file="dge/DGE_ASDvsCTL.xlsx",sheetName = "ASD dge",row.names=TRUE, showNA=FALSE)
xlsx::write.xlsx(sign, file="dge/DGE_ASDvsCTL.xlsx",sheetName = "ASD sign",row.names=TRUE, showNA=FALSE,append=TRUE)
write.table(DGE,"dge/DGE_ASDvsCTL.txt",sep="\t",quote=F)

# Differential expression by region
pd <- demo 
mod <- list()
total <- list()
Y.b <- list()
fitLM <- list()
fitEb <- list()
brDGE <- list()
brDGE_sign <- list()
new_pd <- split(demo,demo$Region) %>% drop.levels()

for (i in 1:11)
        {
            Y.b[[i]]=dat[,colnames(dat) %in% rownames(new_pd[[i]])]
            Y.b[[i]]=Y.b[[i]][,match(rownames(new_pd[[i]]),colnames(Y.b[[i]]))]
            new_pd[[i]]$Region <- NULL 
            new_pd[[i]]$Subject <- NULL
            total[[i]] <- as.formula(paste("~", paste(names(new_pd[[i]]),collapse="+")))
            mod[[i]] <- model.matrix(total[[i]], data = new_pd[[i]])
            fitLM[[i]] = lmFit(Y.b[[i]],mod[[i]],method="robust");
            fitEb[[i]] = eBayes(fitLM[[i]]);
            brDGE[[i]] = topTable(fitEb[[i]], coef = "DiagnosisCTL",number=nrow(Y.b[[i]]),sort.by="none");
            brDGE_sign[[i]] <- brDGE[[i]][brDGE[[i]]$adj.P.Val < 0.05 & abs(brDGE[[i]]$logFC) > 0.3,]
        }

names(brDGE) <- names(new_pd)
names(brDGE_sign) <- names(new_pd)
save(brDGE,brDGE_sign, file="dge/DGE_ASDvsCTL_RegByReg.RData")

for (i in 1:11){
    brDGE_sign[[i]] <- brDGE_sign[[i]] %>%
                            as.data.frame()   %>%
                            rownames_to_column('Gene') 
}

for (i in 1:11){
    brDGE[[i]] <- brDGE[[i]] %>%
                            as.data.frame()   %>%
                            rownames_to_column('Gene') 
}

# Save as excel files
wb <- createWorkbook()
Map(function(data, nameofsheet){     
    addWorksheet(wb, nameofsheet)
    writeData(wb, nameofsheet, data)
}, brDGE_sign, names(brDGE_sign))
saveWorkbook(wb, file = "dge/DGE_ASDvsCTL_RegByReg_Sign.xlsx", overwrite = TRUE)

# Save as excel files
wb <- createWorkbook()
Map(function(data, nameofsheet){     
    addWorksheet(wb, nameofsheet)
    writeData(wb, nameofsheet, data)
}, brDGE, names(brDGE))
saveWorkbook(wb, file = "dge/DGE_ASDvsCTL_RegByReg.xlsx", overwrite = TRUE)

sessionInfo()