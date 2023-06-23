# RNA-seq processing

*This README will be finished soon* 

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
Using featureCouts, we can 
```{bash}
featureCounts -a $GTF -p -B -C -T 15 -o $OUTFILE.csv $BAM1 $BAM2 
```

## 4. Differential Analysis
```{r}
```

