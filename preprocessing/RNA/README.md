# RNA-seq processing

Here, you can find a detailed description of how RNA-seq data was processed in this project.

## Dependencies

* [Trim Galore](https://github.com/FelixKrueger/TrimGalore)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [STAR](https://github.com/alexdobin/STAR/tree/master)
* [sambamba](https://github.com/biod/sambamba)
* [featureCounts](https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html)
* [DESeq2 (R package)](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)


## Summary of the workflow

1. **Trimming and Quality**: Trim Galore, samtools and FastQC
2. **Alignment**: Bowtie2 
3. **Counting**: sambamba and featureCounts
4. **Differential Analysis**: DESeq2


## 1. Trimming and Quality
The first step of the pipeline is the trimming of the reads to remove any possible adapters present in the reads, using trim galore.

```{bash}
trim_galore $FASTQ1 $FASTQ2  --paired --output_dir $FASTQ_DIR --basename $NAME -a $ADAPTERS -A $ADAPTERS --fastqc_args '--outdir $OUTDIR/quality' --cores $THREADS
```

## 2. Alignment
Once the reads are cleaned from adapters we can algin them to the reference genome using STAR.
In this line of code we are converting the SAM file to a BAM file immediately after aligning.
```{bash}
STAR --runThreadN $THREADS --genomeDir $INDICES --readFilesIn $FASTQ1 $FASTQ2 --outFileNamePrefix $OUTNAME --outFilterType BySJout --readFilesCommand zcat --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outStd SAM | samtools view -bS - > $BAMOUT
```
## 3. Counting
We can now quantify the reads falling in each gene of your chosen genome, using its GTF file. 

### a. Sort 
First, we have to sort the BAM file
```{bash}
sambamba sort -t 15 -o /gpfs/projects/bsc08/bsc08471/p53/results/HCT116/Omics/RNA/HCT116_WTNutlin1h/results/RNA/alignment/HCT116_WTNutlin1h_RNA_1.sort.bam /gpfs/projects/bsc08/bsc08471/p53/results/HCT116/Omics/RNA/HCT116_WTNutlin1h/results/RNA/alignment/HCT116_WTNutlin1h_RNA_1.bam
```

### b. Count
Using featureCounts, we can obtain a count matrix of the chosen samples
```{bash}
featureCounts -a $GTF -p -B -C -T 15 -o $OUTFILE.csv $BAM1 $BAM2
```

## 4. Differential Analysis

Finally, we can perform the differential expression analysis between samples of your choice.

### a. Create count matrix
Before performing the differential expression analysis, you have to prepare a count matrix with all the samples and corresponding replicates you will to use.

```{r}
files <- list.files(path = "/path/to/featureCounts/out", pattern = "*.csv$", recursive = T, full.names = T)

countdata <- matrix(data.table::fread(files[1])$Geneid)
colnames(countdata) <- "Geneid"
for (i in files) 
{
  a <- data.table::fread(i)
  a <- a[ ,-(2:6)]
  colnames(a)[-1] <- gsub("\\.[sb]am$", "", basename(colnames(a)[-1]))
  a <- as.matrix(a)
  countdata <- merge(countdata,a, by="Geneid")
}

countdata[,-1] <- apply(countdata[,-1], 2, as.numeric)

data.table::fwrite(countdata, file = "count_table.txt",col.names = T, row.names = F, quote = F, sep = "\t")
```
### b. DEA
Using the complete count matrix, you can now use DESeq2 to perform your analysis
```{r}
# read count matrix
count_matrix <- read.table("path/to/count_matrix.txt",header=T)
colnames(count_matrix) <- c("geneID","Sample1_1", "Sample1_2", "Sample2_1", "Sample2_2")
rownames(count_matrix) <- count_matrix$geneID
count_matrix$geneID <- NULL


annot <- data.frame(sample = colnames(count_matrix))
annot <- cbind(annot, annot %>% separate(col = sample, into = c("condition","replicate"),sep = "_"))
annot$condition <- factor(annot$condition,levels = unique(annot$condition))
annot$condition <- relevel(annot$condition,"Sample1")


dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = annot,
                              design = ~ condition)

dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

res <- results(dds,contrast=c("condition","Sample1","Sample2"))

```


