library(tidyverse)

asd <- read.table("ASD_fMRI_DevClustered.txt",header=T)
krien <- read.table("Krienen.txt",header=T)

df <- merge(asd,krien,by.x="Gene",by.y="Genes",all=F) %>%
  	  pivot_longer(cols = c("PVALB","SST","ID2","VIP"), names_to = "Cell", values_to = "Exp") %>%
  	  mutate(logExp = log2(Exp + 1),Species = as.factor(Species),Cell = as.factor(Cell),Direction = as.factor(Direction)) %>%
  	  mutate(Species = fct_relevel(Species, "Human", "Macaque","Marmoset","ferret"))


my_comparisons <- list( c("Human", "Macaque"), c("Human", "Marmoset"), c("Human", "ferret") )

pdf("Adult_Krienen.pdf",width=3,height=5)
ggboxplot(df %>% filter(Direction == "Adult"), "Species", "logExp", color = "Species",
 palette = "jco",outlier.shape = NA) + 
 facet_wrap(.~Cell,ncol=2,nrow=2) + 
  stat_compare_means(comparisons = my_comparisons,method = "t.test",
                     hide.ns = TRUE) +
theme_classic() +
  xlab("") + 
  rotate_x_text(angle = 45) + 
  ylab("Expression")+
  ylim(0,12) +
 theme(legend.position="none")
 dev.off()


pdf("EarlyDev_Krienen.pdf",width=3,height=5)
ggboxplot(df %>% filter(Direction == "EarlyDev"), "Species", "logExp", color = "Species",
 palette = "jco",outlier.shape = NA) + 
 facet_wrap(.~Cell,ncol=2,nrow=2) + 
  stat_compare_means(comparisons = my_comparisons,method = "t.test",
                     hide.ns = TRUE)  +
theme_classic() +
  xlab("") + 
  rotate_x_text(angle = 45) + 
  ylab("Expression")+
  ylim(0,12) +
 theme(legend.position="none")
 dev.off()


pdf("Stable_Krienen.pdf",width=3,height=5)
ggboxplot(df %>% filter(Direction == "Stable"), "Species", "logExp", color = "Species",
 palette = "jco",outlier.shape = NA) + 
 facet_wrap(.~Cell,ncol=2,nrow=2) + 
  stat_compare_means(comparisons = my_comparisons,method = "t.test",
                     hide.ns = TRUE) +
theme_classic() +
  xlab("") + 
  rotate_x_text(angle = 45) + 
  ylab("Expression")+
  ylim(0,12) +
 theme(legend.position="none")
 dev.off()