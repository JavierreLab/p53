# ChIP-seq processing

Here, you can find a detailed description of how ChIP-seq data was processed in this project.

## Dependencies

* [Trim Galore](https://github.com/FelixKrueger/TrimGalore)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [samtools](http://www.htslib.org/download/)
* [sambamba](https://github.com/biod/sambamba)
* [deepTools](https://deeptools.readthedocs.io/en/latest/)
* [macs2](https://github.com/macs3-project/MACS)
* [csaw (R package)](https://bioconductor.org/packages/release/bioc/html/csaw.html)
* [GenomicRanges (R package)](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
* [DESeq2 (R package)](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)


## Summary of the workflow

1. **Trimming and Quality**: Trim Galore and FastQC
2. **Alignment**: Bowtie2
3. **Filtering**: samtools and sambamba
4. **Coverage**: deepTools
5. **Peak Calling**: Macs2
6. **Generating Peakmatrix**: csaw and GenomicRanges
7. **Generating Background Matrix**: csaw and GenomicRanges
8. **Differential Analysis**: DESeq2


A detailed description of how steps 1-4 were performed can be found [here](https://github.com/JavierreLab/liCHiC/tree/main/3.ChIPseq%20Processing), the github page associated to Tomás-Daza, L *et al.* Low input capture Hi-C (liCHi-C) identifies promoter-enhancer interactions at high-resolution. *Nature Communications* **14**, 268 (2023). [https://doi.org/10.1038/s41467-023-35911-8](https://doi.org/10.1038/s41467-023-35911-8)

## 5. Peak Calling

This step computes significants peaks detected in the experiment. Different variations used in this project are shown below.
```{bash}
# Calling peaks in default (narrow) mode for histone mark H3K27ac.
macs2 callpeak -t $BAM -f $FORMAT -c $BAMINPUT -g $ORG -B -n ${PREFIX} --outdir $OUTDIR/peaks

# Calling peaks in broad mode for histone mark H3K27ac.
macs2 callpeak --broad -t $BAM -f $FORMAT -c $BAMINPUT -g $ORG -B -n ${PREFIX} --outdir $OUTDIR/peaks

# Calling merge/consensus peaks for one sample with many replicates
macs2 callpeak -f $FORMAT -t $BAM1 $BAM2 $BAM3  -c $BAMINPUT1 $BAMINPUT2 $BAMINPUT3  -g $ORG -B -n ${PREFIX} --outdir $OUTDIR/peaks
```

## 6. Generating Peakmatrix
This step allows us to generate a peakmatrix of all our samples. By doing this, we can quantify and compare peaks between samples.

### a. Choose samples of interest
For this, you need to have a complete vector of paths to your sample's bam files (obtained afters step 3) and corresponding peak files (obtained after step 5)
```{r}
library(csaw)
library(GenomicRanges)
library(tidyverse)

path_bams <- arg[1] # list of paths to your bam files of interest
path_peaks <- arg[2] # list of paths to your corresponding peak files of interest
```
### b. Read in peaks
Here, you read all peak files and generate one object with the complete list of all significant peaks found in all samples of interest
```{r}
# read in peak files, covert to GRanges objects
peaks <- GRanges()
for (i in peak_files) 
{
  a <- data.table::fread(i)[,1:4]
  aGR <- makeGRangesFromDataFrame(a, seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = T)
  peaks <- c(peaks, aGR)
}
colnames(elementMetadata(peaks)) <- "name"
peaks$ID <- 1:length(peaks)
```

### c. Disjoin peaks
The next step is to perform a disjoin of the peaks. This means that two overlapping peaks are split into three regions, with region 1 corresponding to the region unique to peak1, regions 2 being the region where both peaks overlap, and region 2 correspoing to the region unique to peak3. The figure below illustrates this.  

![Disjoin](disjoin.png)

```{r}
# disjoin peaks
peak_disj <- disjoin(peaks,with.revmap=T)
```

### d. Count reads in regions
In this step we can quantify the number of reads falling in each region per bam file (per sample)
```{r}
# count reads falling in each region of peak_disj for each sample in bam_files
param <- readParam()

# counting
peak_count <- regionCounts(bam.files = bam_files, regions = peak_disj,
                           ext = NA, param = param, 
                           BPPARAM = BatchtoolsParam(workers = length(bam_files)*2))   
```

### e. Format and save peakmatrix
Finally, we can format the resulting quantification and save as a peakmatrix, where each row reflects a region and each column represents a bam file (sample).
```{r}
peak_df <- as.data.frame(assay(peak_count))
colnames(peak_df) <- colData(peak_count)$bam.files
peak_df$region <- paste(seqnames(rowRanges(peak_count)),ranges(rowRanges(peak_count)), sep = ":")

data.table::fwrite(peak_df[,c(ncol(peak_df),1:(ncol(peak_df)-1))], file=paste0(out_dir,"/peak_matrix_counts.tsv"), row.names = F, col.names = T, quote = F, sep = "\t")
```

## 7. Generating Background Matrix
Repeat steps 6a-c first

### d. Count background reads
Quantify background noise in 10kb bins
```{r}

# discard disjoined peaks from background counting
param <- readParam(discard=peak_disj)

# count background noise
bg_count <- windowCounts(bam.files = bam_files, width = 10000, spacing = 10000,
                         ext = NA, filter = 0, param = param, 
                         BPPARAM = BatchtoolsParam(workers = length(bam_files)*2))
```
### e. Format and save background matrix
Finally, we can format the resulting bacground quantification and save as a background matrix, where each row reflects a 10kb bin and each column represents a bam file (sample).
```{r} 
bg_df <- as.data.frame(assay(bg_count))
colnames(bg_df) <- colData(bg_count)$bam.files
bg_df$region <- paste(seqnames(rowRanges(bg_count)),ranges(rowRanges(bg_count)), sep = ":")
 
data.table::fwrite(bg_df[,c(ncol(bg_df),1:(ncol(bg_df)-1))], file=paste0(out_dir,"/background_counts.tsv"), row.names = F, col.names = T, quote = F, sep = "\t")
```

## 8. Differential Analysis
In this final step, we perform the differential analysis of a histone mark. The output will be a table with regions found with a significant differential deposition of a histone mark (gain or loss). 

### a. Estimate size factors of background matrix  
Size factors (normalization factors) of bacground matrix have to be computed in order to perform the background correction to our peak matrix for proper normalisation before differential analysis.
```{r}
library(DESeq2)

# read background matrix
background_matrix <- read.table("path/to/background_matrix.txt",header=T)
colnames(background_matrix) <- c("regions","Sample1_1", "Sample1_2", "Sample2_1", "Sample2_2")
rownames(background_matrix) <- background_matrix$regions
background_matrix$regions <- NULL

# prepare metadata
coldata <- data.frame(id=c("Sample1_1", "Sample1_2", "Sample2_1", "Sample2_2") ,
                      tissue=c("Sample1", "Sample1", "Sample2", "Sample2"),
                      rep=c(1, 2, 1, 2))
rownames(coldata) <- coldata$id


dds <- DESeqDataSetFromMatrix(countData = background_matrix,
                              colData = coldata,
                              design = ~ tissue)

# estimate size factors
dds <- estimateSizeFactors(dds)
sf <- sizeFactors(dds)
```
### b. Background correction
Correct peakmatrix for background signal using the background matrix's size factors.
```{r}

# read peakmatrix
count_matrix <- read.table("path/to/peak_matrix.txt",header=T)
colnames(count_matrix) <- c("regions","Sample1_1", "Sample1_2", "Sample2_1", "Sample2_2")
rownames(count_matrix) <- count_matrix$regions
count_matrix$regions <- NULL

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ tissue)

# apply pre-compute size factors with background regions - background correction
sizeFactors(dds) <- sf

```
### c. Differential Analysis  
Finally, we can perform the differential analysis between samples. Remember this has to be done with more than one replicate per sample, so if you want to perform this step, you have to generate the peak and background matrices with this in mind.
```{r}
dds <- DESeq(dds)
res <- results(dds,contrast=c("tissue","Sample1","Sample2"))
```






