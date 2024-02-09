<img src="https://github.com/JavierreLab/p53/assets/86778675/fca5a6f5-e09a-44a1-a12d-9e1fcc0d6299"  width="300" ALIGN="right">

# p53 rapidly restructures 3D chromatin organization to trigger a transcriptional response
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8075024.svg)](https://doi.org/10.5281/zenodo.8075024)

This repository contains the necessary data and scripts to reproduce the main figures of the  manuscript.

The data needed to execute the preprocessing step is available at [GEO](https://www.ncbi.nlm.nih.gov/geo/) under the accession GSMxxxxx. In the corresponding `preprocessing` folder are listed and organized by omics the different scripts used to process our data from the trimming and mapping of FASTQ files to the call of significant biological signals (peaks, interactions...).

Nonetheless, the `data` folder allows to skip the preprocessing step as it contains the supplementary data tables and essential processed data needed to generate most figures presented the manuscript.

The `figures` folder contains script and folder to (re)generate stats and figures from the manuscript.  

## Dependencies

Specific dependencies are listed in each `preprocessing` section (see below), and in the `figures` section.

## Preprocessing

| Datatype | Short description |
|----------|------|
| [ChIP-Seq](https://github.com/JavierreLab/p53/tree/main/preprocessing/ChIP) |  Trimming, mapping, peak-calling, and differential analysis  |
| [Hi-C](https://github.com/JavierreLab/p53/tree/main/preprocessing/HiC) |  Mapping, filtering, normalization, compartment- and TAD-calling  |
| [PCHi-C](https://github.com/JavierreLab/p53/tree/main/preprocessing/PCHiC)    |  Mapping and significant interactions calling |
| [RNA-seq](https://github.com/JavierreLab/p53/tree/main/preprocessing/RNA)  |  Trimming, mapping and differential expression analysis  |
