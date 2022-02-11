# Enrichment.r
# -g = Two columns table of input genes with specific association from your study
# -l = list of two columns tables with gene - disease association. E.g. Gene1 SYN
# -p = make a bubble chart with OR and -log10(FDR)
# -b = background (protein coding = 19776, brain expressed = 15585)
# -o = output label for statistics and viz
# -W/-H = width/height of the plot. 

mkdir enrichments_dcgenes/

cp rawdata/geneset_for_enrichment/*.RData enrichments_dcgenes/
cp UTILS/Enrichment.r enrichments_dcgenes/
cp integrative_results/ASD_fMRI_Genes.txt enrichments_dcgenes/
cp integrative_results/ASD_fMRI_Genes2.txt enrichments_dcgenes/

cd enrichments_dcgenes/

mkdir STATS/

# Allen scRNA
Rscript Enrichment.r -g ASD_fMRI_Genes.txt -l Allen_MultiReg_CellMarkers_GeneSet.RData -p -b 15585 -o STATS/Allen_scRNA -W 5 -H 5
Rscript Enrichment.r -g ASD_fMRI_Genes.txt -l Allen_MultiReg_FULL_CellMarkers_GeneSet.RData -p -b 15585 -o STATS/Allen_scRNA_Full -W 6 -H 15
Rscript Enrichment.r -g ASD_fMRI_Genes.txt -l GeneSets_BBLake_2018.RData -p -b 15585 -o STATS/BBlake_scRNA -W 3 -H 6
Rscript Enrichment.r -g ASD_fMRI_Genes.txt -l Allen_MultiReg_LayersMarkers_GeneSet.RData -p -b 15585 -o STATS/Allen_scRNA_Full_Layers -W 4 -H 10

# Neuropsy
Rscript Enrichment.r -g ASD_fMRI_Genes.txt -l ASD_SFARI.RData -p -b 15585 -o STATS/DISORDERS -W 5 -H 5
Rscript Enrichment.r -g ASD_fMRI_Genes.txt -l PsychENCODE_DEGs.RData -p -b 15585 -o STATS/PSY_DEGS -W 5 -H 5
Rscript Enrichment.r -g ASD_fMRI_Genes.txt -l PsychEncode_Modules.RData -p -b 15585 -o STATS/PSY_MODS -W 3 -H 5
Rscript Enrichment.r -g ASD_fMRI_Genes.txt -l ASD_SingleCell_DEGs.RData -p -b 15585 -o STATS/ASD_scRNA -W 5 -H 5
Rscript Enrichment.r -g ASD_fMRI_Genes.txt -l CMC_DEGs.RData -p -b 15585 -o STATS/CMC_DEGS -W 5 -H 5

Rscript Enrichment.r -g ASD_fMRI_Genes.txt -l GeneSets_CHD8_Targets.RData -p -b 15585 -o STATS/CHD8_Targets -W 5 -H 5



# Develop
Rscript Enrichment.r -g ASD_fMRI_Genes.txt -l GeneSets_Parikshak.RData -p -b 15585 -o STATS/DEVBRAIN_Parikshak -W 5 -H 5
Rscript Enrichment.r -g ASD_fMRI_Genes.txt -l GeneSets_Polioudakis2020.RData -p -b 15585 -o STATS/DEVBRAIN_Polioudakis -W 5 -H 5

# DEGS by Reg
Rscript Enrichment.r -g ASD_fMRI_Genes.txt -l DEGs_ByReg.RData -p -b 15585 -o STATS/REGIONAL_DGE -W 5 -H 5

# DEGS by Reg
Rscript Enrichment.r -g ASD_fMRI_Genes.txt -l fMRI_Genes.RData -p -b 15585 -o STATS/POS_NEG_TopRsq -W 5 -H 5

Rscript Enrichment.r -g ASD_fMRI_Genes.txt -l GeneSets_Modules_Jillian.RData -p -b 15585 -o STATS/Modules_Jillian -W 10 -H 5

Rscript Enrichment.r -g ASD_fMRI_Genes2.txt -l GeneSets_fMRI.RData -p -b 15585 -o STATS/Previous_fMRI -W 5 -H 5

rm *.RData