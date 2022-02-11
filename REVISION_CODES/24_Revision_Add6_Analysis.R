library(Hmisc)
library(tidyverse)
library(ggpubr)
library(here)

load(here("integrative_results","DiffCor_Significant.RData"))

Neuron <- read.table("rawdata/fALFF_genes_Pearson.txt",header=T)
Neuron <- data.frame(Gene=Neuron$fALFF_genes, Class=rep("Neuron",nrow(Neuron)))
Rich <- read.table("rawdata/Richiardi_science.txt",header=T)

universe <- read.table("UTILS/BrainExpressed_Background.txt",header=T)

P=1000  ## select number of permutation resamples

Y.a <- replicate(P,sample(universe$x,nrow(Rich),replace = TRUE),simplify=FALSE)
Y.b <- replicate(P,sample(universe$x,nrow(Neuron),replace = TRUE),simplify=FALSE)

# Delta Simulation
boot.a <- list()
for(i in 1:P) {
    boot.a[[i]] <- length(intersect(DiffCor_Sign_Final$Gene,Y.a[[i]]))
}

boot.b <- list()
for(i in 1:P) {
    boot.b[[i]] <- length(intersect(DiffCor_Sign_Final$Gene,Y.b[[i]]))
}


pdf("revision/Permutation_DiffCor_Rich.pdf",width=4,height=4)
intercept <- length(intersect(DiffCor_Sign_Final$Gene,Rich$Gene))
df.a <- do.call(rbind,boot.a) %>% as.data.frame()
k <- sum(df.a[,1] >= intercept)
p <- round(zapsmall(binconf(k, length(df.a[,1]), method='exact'))[3],3)
gghistogram(df.a, 
      x = "V1", 
      color = "grey60",
      add = "mean", rug = FALSE,
      bins = 15)+
      geom_vline(xintercept=intercept, colour="black",size=0.5) +
      annotate(geom="text", x=intercept+1, y=300, label=paste0("p: ",p)) +
      xlim(0,15) + ylim(0,350) +
      xlab("# of Genes")
dev.off()

pdf("revision/Permutation_DiffCor_Wang.pdf",width=4,height=4)
intercept <- length(intersect(DiffCor_Sign_Final$Gene,Neuron$Gene))
df.b <- do.call(rbind,boot.b) %>% as.data.frame()
k <- sum(df.b[,1] >= intercept)
p <- round(zapsmall(binconf(k, length(df.b[,1]), method='exact'))[3],3)
gghistogram(df.b, 
      x = "V1", 
      color = "purple",
      add = "mean", rug = FALSE,
      binwidth = 0.5)+
      geom_vline(xintercept=intercept, colour="black",size=0.5) +
      annotate(geom="text", x=intercept+1, y=300, label=paste0("p: ",p)) +
      xlim(0,8) + ylim(0,500) +
      xlab("# of Genes")
dev.off()


# CTL 

# Delta Simulation
boot.a <- list()
for(i in 1:P) {
    boot.a[[i]] <- length(intersect(CTL_Sign_Final$Gene,Y.a[[i]]))
}

boot.b <- list()
for(i in 1:P) {
    boot.b[[i]] <- length(intersect(CTL_Sign_Final$Gene,Y.b[[i]]))
}


pdf("revision/Permutation_CTLsign_Rich.pdf",width=4,height=4)
intercept <- length(intersect(CTL_Sign_Final$Gene,Rich$Gene))
df.a <- do.call(rbind,boot.a) %>% as.data.frame()
k <- sum(df.a[,1] >= intercept)
p <- round(zapsmall(binconf(k, length(df.a[,1]), method='exact'))[3],3)
gghistogram(df.a, 
      x = "V1", 
      color = "grey60",
      add = "mean", rug = FALSE)+
      geom_vline(xintercept=intercept, colour="black",size=0.5) +
      annotate(geom="text", x=intercept+1, y=300, label=paste0("p: ",p)) +
      xlim(0,70) + ylim(0,350) +
      xlab("# of Genes")
dev.off()

pdf("revision/Permutation_CTLsign_Wang.pdf",width=4,height=4)
intercept <- length(intersect(CTL_Sign_Final$Gene,Neuron$Gene))
df.b <- do.call(rbind,boot.b) %>% as.data.frame()
k <- sum(df.b[,1] >= intercept)
p <- round(zapsmall(binconf(k, length(df.b[,1]), method='exact'))[3],3)
gghistogram(df.b, 
      x = "V1", 
      color = "purple",
      add = "mean", rug = FALSE)+
      geom_vline(xintercept=intercept, colour="black",size=0.5) +
      annotate(geom="text", x=intercept+1, y=300, label=paste0("p: ",p)) +
      xlim(0,40) + ylim(0,350) +
      xlab("# of Genes")
dev.off()


