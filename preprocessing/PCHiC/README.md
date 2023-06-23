# PCHi-C processing

Here, you can find a detailed description of how PCHi-C data was processed in this project.

## Dependencies

* [HiCUP](https://www.bioinformatics.babraham.ac.uk/projects/hicup/)
* [Chicago (R package)](https://www.bioconductor.org/packages/release/bioc/html/Chicago.html)
* [bowtie2](https://github.com/BenLangmead/bowtie2)


## Summary of the workflow

1. **Mapping and Filtering**: HiCUP
2. **Capture Efficiency**: HiCUP miscellaneous script
3. **Interaction Calling**: Chicago
4. **Weight Recalibration of Interaction Calling**: Chicago



A detailed description of how steps 1-3 were performed can be found [here](https://github.com/JavierreLab/liCHiC/tree/main/1.liCHiC%20Processing), the github page associated to Tomás-Daza, L *et al.* Low input capture Hi-C (liCHi-C) identifies promoter-enhancer interactions at high-resolution. *Nature Communications* **14**, 268 (2023). [https://doi.org/10.1038/s41467-023-35911-8](https://doi.org/10.1038/s41467-023-35911-8)

## 4. Weight Recalibration of Interacting Calling

Default parameters for interaction calling are calibrated on high-confidence calls from seven human Macrophage data sets (i.e. interactions that pass our p-value threshold in all seven samples). Provided that your cell type is not too dissimilar to these calibration data, it should be fine to leave the parameters at their default settings. However, if your data set is from an unusual cell type, you may wish to recalibrate these parameters using data from cell types similar to yours. 
To do this, you should used the *fitDistCurve* function from ChicagoTools. Further explanation on how this can be done is found in [Chicago's bitbucket](https://bitbucket.org/chicagoTeam/chicago/src/master/chicagoTools/)

In our case, our datasets are from HCT116 cell line, a cancer cell line. For this reason, we decided to recalibrate the parameters used to our control data. 

