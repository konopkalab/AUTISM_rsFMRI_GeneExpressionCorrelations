Association between resting-state functional brain connectivity and gene expression is altered in autism spectrum disorder 
==========================

This repository contains analysis code for the rs-fMRI ~ gene expression project carried out by researchers at the [Konopka Lab, UTSW](http://konopkalab.org/) and [Montillo Lab, UTSW](https://aamontillo.net/)

## Cite this

If you use anything in this repository please cite the following publication:

Pre-print URL: https://www.medrxiv.org/content/10.1101/2021.01.07.21249281v1


## Files

| directory | contents | code |
| --------- | -------- | -------- |
| [`processing_qc`](processing_qc/) | Output data from initial processing and quality check. | 01_Data_processing_QC.R |
| [`DGE`](DGE/) | Output data from DGE analysis. | 02_DGE_Analysis.R |
| [`MAGMA`](MAGMA/) | Output data from MAGMA gene set analysis. | 03_MAGMA.sh |
| [`Shiny_App`](Shiny_App/) | Input data to visualize the data. | HippoAxisSeq_App.R |

## Explore the data

We have provided an interactive web app that allow you to explore the data at single nucleus level. 

* https://human-hippo-axis.cells.ucsc.edu/

* Shiny app (Please download the data and build it using the HippoAxisSeq_App.R code).

![](HippoAxisSeq.gif)]
