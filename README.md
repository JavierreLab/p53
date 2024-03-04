
# p53 rapidly restructures 3D chromatin organization to trigger a transcriptional response
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8075024.svg)](https://doi.org/10.5281/zenodo.8075024)

<img src="https://github.com/JavierreLab/p53/assets/86778675/fca5a6f5-e09a-44a1-a12d-9e1fcc0d6299"  width="300" ALIGN="right">

This repository contains the necessary data and scripts to reproduce the main figures of the  manuscript.

The data needed to execute the preprocessing step is available at [GEO](https://www.ncbi.nlm.nih.gov/geo/) under the accession [GSE235947](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235947). In the corresponding `preprocessing` folder are listed and organized by omics the different scripts used to process our data from the trimming and mapping of FASTQ files to the call of significant biological signals (peaks, interactions...).

Nonetheless, the `data` folder allows to skip the preprocessing step as it contains the supplementary data tables and essential processed data needed to generate most figures presented the manuscript.

The `figures` folder contains script and folder to (re)generate stats and figures from the manuscript.  

## Dependencies

Specific dependencies are listed in each [preprocessing](https://github.com/JavierreLab/p53/tree/main/preprocessing) section:

| Datatype | Short description |
|----------|------|
| [ChIP-Seq](https://github.com/JavierreLab/p53/tree/main/preprocessing/ChIP) |  Trimming, mapping, peak-calling, and differential analysis  |
| [Hi-C](https://github.com/JavierreLab/p53/tree/main/preprocessing/HiC) |  Mapping, filtering, normalization, compartment- and TAD-calling  |
| [PCHi-C](https://github.com/JavierreLab/p53/tree/main/preprocessing/PCHiC)    |  Mapping and significant interactions calling |
| [RNA-seq](https://github.com/JavierreLab/p53/tree/main/preprocessing/RNA)  |  Trimming, mapping and differential expression analysis  |

and in the [figures](https://github.com/JavierreLab/p53/tree/main/figures/) section:

|||
|---|---|
| [Figure 1](https://github.com/JavierreLab/p53/tree/main/figures/Figure_1) |  A/B compartment changes along p53 activation. |
| [Figure 2](https://github.com/JavierreLab/p53/tree/main/figures/Figure_2) |  TAD dynamics along p53 activation. |
| [Figure 3](https://github.com/JavierreLab/p53/tree/main/figures/Figure_3) |  p53 binding to chromatin leads to direct changes in 3D genome topology. |
| [Figure 4](https://github.com/JavierreLab/p53/tree/main/figures/Figure_4) |  Promoter interactions shift during p53 activation. |
| [Figure 5](https://github.com/JavierreLab/p53/tree/main/figures/Figure_5) |  Identification of new p53 distal target genes through enhancer-promoter interactions. |
| [Figure 6](https://github.com/JavierreLab/p53/tree/main/figures/Figure_6) |  p53 drives the formation of promoter-enhancer interactions or uses pre-established ones to control gene transcription over distance. |
| [Figure 7](https://github.com/JavierreLab/p53/tree/main/figures/Figure_7) |  Cohesin depletion impedes p53-mediates transcriptional response.  |
| [Supplementary](https://github.com/JavierreLab/p53/tree/main/figures/Supplementary_Figures_bulk) |  Scripts related to supplementary panels. |
