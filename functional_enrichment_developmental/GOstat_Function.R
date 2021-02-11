
dogo <- function(names,universe,species="human", goP = 0.01, 
	cond=FALSE, ontology = "BP"){
    if(species=="human"){
		golib="org.Hs.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Hs.egSYMBOL2EG
  } else  if (species == "mouse") {
		golib="org.Mm.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Mm.egSYMBOL2EG
	}
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA))
  x=x[!is.na(x)]
  Universe=unlist(mget(as.character(universe),gomap,ifnotfound = NA))
 Universe=unique(c(Universe[!is.na(Universe)],unique(x)))

  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = ontology, pvalueCutoff = goP, conditional = cond,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)

}

# How to run
files=list.files(pattern="ASD_fMRI")
mod=as.data.frame(lapply(files,read.table,sep="\t",header=T)[[1]])
colnames(mod)[2] <- "Class"
uni <- read.table("universe.txt",header=F)
uni <- as.character(uni$V1)

list=list()
df=list()
for (net in levels(mod$Class)) {
list[[net]]=mod[mod$Class == net,]
df[[net]]=dogo(list[[net]]$Gene,uni,species="human", goP = 1,cond=FALSE, ontology = "BP")
df[[net]]$adj=p.adjust(df[[net]]$Pvalue,"BH")
df[[net]]=df[[net]][df[[net]]$Pvalue < 0.05 & df[[net]]$Count > 5,]
}
dir.create("BP/")
net=levels(mod$Class)
for (i in 1:length(df)) write.table(df[[i]], file = paste("BP/BP_", net[[i]], ".txt",sep = ""),sep="\t",quote=F)


list=list()
df=list()
for (net in levels(mod$Class)) {
list[[net]]=mod[mod$Class == net,]
df[[net]]=dogo(list[[net]]$Gene,uni,species="human", goP = 1,cond=FALSE, ontology = "MF")
df[[net]]$adj=p.adjust(df[[net]]$Pvalue,"BH")
df[[net]]=df[[net]][df[[net]]$Pvalue < 0.05 & df[[net]]$Count > 5,]
}
dir.create("MF/")
net=levels(mod$Class)
for (i in 1:length(df)) write.table(df[[i]], file = paste("MF/MF_", net[[i]], ".txt",sep = ""),sep="\t",quote=F)


list=list()
df=list()
for (net in levels(mod$Class)) {
list[[net]]=mod[mod$Class == net,]
df[[net]]=dogo(list[[net]]$Gene,uni,species="human", goP = 1,cond=FALSE, ontology = "CC")
df[[net]]$adj=p.adjust(df[[net]]$Pvalue,"BH")
df[[net]]=df[[net]][df[[net]]$Pvalue < 0.05 & df[[net]]$Count > 5,]
}
dir.create("CC/")
net=levels(mod$Class)
for (i in 1:length(df)) write.table(df[[i]], file = paste("CC/CC_", net[[i]], ".txt",sep = ""),sep="\t",quote=F)
