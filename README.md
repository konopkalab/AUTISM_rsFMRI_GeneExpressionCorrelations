Association between resting-state functional brain connectivity and gene expression is altered in autism spectrum disorder 
==========================

This repository contains analysis code for the rs-fMRI ~ gene expression project carried out by researchers at the [Konopka Lab, UTSW](http://konopkalab.org/) and [Montillo Lab, UTSW](https://aamontillo.net/)

## Cite this

If you use anything in this repository please cite the following publication:

Pre-print URL: https://www.medrxiv.org/content/10.1101/2021.01.07.21249281v1


## Files

| directory | contents | code |
| --------- | -------- | -------- |
| [`rawdata`](rawdata/) | Input/Output data of the initial processing and quality check. | 01_Data_processing_QC.R \ 02_Data_Preprocess_QC_Visualization.R \ 03_Imaging_Processing.R|
| [`dge`](dge/) | Output of the Differential expression analysis. | 01_Data_Preprocess_QC.R|
| [`Permuted_Matches_Outputs`](Permuted_Matches_Outputs/) | Output data of the correlative analysis. | 04_Permuted_Matches_fALFF.R \ 05_fALFF_LOR_Analysis.R \ 06_Permuted_Matches_ReHo.R \ 07_ReHo_LOR_Analysis.R|
| [`imaging_visualizations`](imaging_visualizations/) | Output data of the visualizations for rs-fMRI data. | 08_Imaging_Visualizations.R |
| [`integrative_results`](integrative_results/) | Output of the integration analysis. | 09_Data_Integration.R |
| [`integrative_visualizations`](integrative_visualizations/) | Output of the visualizations for the integrated data and leave one region out analysis. | 10_Data_Visualization.R \ 12_Leave-one-region-out_Contribution.R|
| [`motif_enrichment`](motif_enrichment/) | Output of the motif enrichment analysis. | 11_Motif_Enrichment.R |
| [`enrichments_developmental`](enrichments_developmental/) | Output of the gene set enrichment analysis for developmental clusters. | 13_Run_Enrichment.sh |
| [`enrichments_dcgenes`](enrichments_dcgenes/) | Output of the gene set enrichment analysis for the differentially correlated genes. | 14_Run_Enrichment_DCgenes.sh |
| [`supp_tables`](supp_tables/) | Output of the database creation. | 15_Database.R |
| [`networking`](networking/) | Output of the PPI network analysis. | 17_PPI_Networks.R|
