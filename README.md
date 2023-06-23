TODO: make Zenodo (only when public)

<img src="https://github.com/JavierreLab/p53/assets/106498/c485c9ec-c9d2-4838-89a5-89c783b1ac26" width="250" ALIGN="left">

# p53 rapidly restructures 3D chromatin organization to trigger a transcriptional response

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
