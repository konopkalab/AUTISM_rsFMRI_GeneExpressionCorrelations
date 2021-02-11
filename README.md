Association between resting-state functional brain connectivity and gene expression is altered in autism spectrum disorder 
==========================

This repository contains analysis code for the rs-fMRI ~ gene expression project carried out by researchers at the [Konopka Lab, UTSW](http://konopkalab.org/) and [Montillo Lab, UTSW](https://aamontillo.net/)

## Cite this

If you use anything in this repository please cite the following publication:

Pre-print URL: https://www.medrxiv.org/content/10.1101/2021.01.07.21249281v1


## Files

| directory | contents | code |
| --------- | -------- | -------- |
| [`rawdata`](rawdata/) | Input/Output data from initial processing and quality check. | 01_Data_processing_QC.R \ 02_Data_Preprocess_QC_Visualization.R \ 03_Imaging_Processing.R|
| [`Permuted_Matches_Outputs`](Permuted_Matches_Outputs/) | Output data from correlative analysis. | 04_Permuted_Matches_fALFF.R \ 05_fALFF_LOR_Analysis.R \ 06_Permuted_Matches_ReHo.R \ 07_ReHo_LOR_Analysis.R|
| [`MAGMA`](MAGMA/) | Output data from MAGMA gene set analysis. | 03_MAGMA.sh |
| [`Shiny_App`](Shiny_App/) | Input data to visualize the data. | HippoAxisSeq_App.R |
