<img src="https://github.com/JavierreLab/p53/assets/86778675/fca5a6f5-e09a-44a1-a12d-9e1fcc0d6299"  width="300" ALIGN="right">

# p53 rapidly restructures 3D chromatin organization to trigger a transcriptional response
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8075024.svg)](https://doi.org/10.5281/zenodo.8075024)

This repository contains the necessary data and scripts to reproduce the main figures of the  manuscript.

The data needed to execute the preprocessing step is available at [GEO](https://www.ncbi.nlm.nih.gov/geo/) under the accession GSMxxxxx. In the corresponding `preprocessing` folder are listed and organized by omics the different scripts used to process our data from the trimming and mapping of FASTQ files to the call of significant biological signals (peaks, interactions...).

Nonetheless, the `data` folder allows to skip the preprocessing step as it contains the supplementary data tables and essential processed data needed to generate most figures presented the manuscript.

The `figures` folder contains script and folder to (re)generate stats and figures from the manuscript.  

## Prerequisites

### python libraries

 - scipy
 - numpy
 - mapltolib
 - TADbit (https://github.com/fransua/tadbit/tree/p53_javierre)
 - bioframe

### R libraries

## Preprocessing

| Datatype | link |
|----------|------|
| ChIP-Seq | |
| Hi-C     | |
| PCHiC    | |
| RNA-seq  | |
