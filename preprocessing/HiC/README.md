# Hi-C processing

*This README will be finished soon* 

Here, you can find a detailed description of how RNA-seq data was processed in this project.

## Dependencies

* [TADbit](https://github.com/fransua/tadbit/tree/p53_javierre) ([Serra et al. 2017](https://doi.org/10.1371/journal.pcbi.1005665))
* [GEM3](https://github.com/smarco/gem3-mapper) ([Marco-Sola et al. 2012](https://doi.org/10.1038/nmeth.2221))

## Summary of TADbit workflow

1. **Alignment**: "restriction-enzyme aware" alignment with GEM3 
2. **Parsing**: pull sequencing lanes and replicates together
3. **Filtering**: remove non-canonical di-tags
4. **Normalizing**: remove Hi-C specific biases
5. **Compartments calling**: call for A/B compartments and manual curation
6. **TAD calling**: call for Topologically Associating Domains (TADs) and insulation score


## 1. Alignment

For each sample TADbit was used to map with [GEM3](https://github.com/smarco/gem3-mapper) the trimmed FASTQ files. 

Here an example with the Nutlin1a-1h sample (one command per read-end):

```
tadbit map -w HCT116_WT_NUT_1h --mapper gem --mapper_binary gem3-mapper --fastq HCT116_HiC_Nut1h_BR1_1.fq.gz --index hg19.gem3 --renz HindIII --read 1 -C 24 --noX  --mapper_param "--alignment-local-min-identity 15"
tadbit map -w Cabrera_2022_HCT116_WT_DMSO --mapper gem --mapper_binary gem3-mapper --fastq HCT116_HiC_Nut1h_BR1_1.fq.gz --index hg19.gem3 --renz HindIII --read 1 -C 24 --noX  --mapper_param "--alignment-local-min-identity 15"
```
*Note: TADbit uses fragment based mapping by default. Full length read end will be mapped first. Then, for unmapped reads, TADbit will search for restriction-ezmyme specific ligation site and map the corresponding read-end fragments independently.*

## 2. Parsing

This step will pull together all mapped results from each read end for the different mapping steps, and for different biological replicates (if applicable).

The command following Nutlin1a-1h example would be:

```
tadbit parse -w HCT116_WT_NUT_1h --genome hg19.fa --noX --compress_input
```

## 3. Filtering

This step consists of two substeps:
  - Finding the intersection between mapped reads in both ends. Forming a di-tags for each read mapped in both ends (*Note that TADbit allows the generation of multiple contacts, not only pairwise, as fragment read-ends by ligation site may result in mapping more than 2 fragments per read-end pair*).
  - Filtering di-tags accoring to a list of parameters:
    1. self-circle (both read-ends in the same restriction fragment mapped in opposite strands facing outwards)
    2. dangling-end (both read-ends in the same restriction fragment mapped in opposite strands facing inwards)
    3. error (both read-ends in the same restriction fragment mapped in same strands)
    4. extra dangling-end (read-ends in different restriction fragment but very close* and mapped in opposite strands facing inwards)
    5. too close from RES (not applied here)
    6. too short (one of the read-end mapped in a restriction-enzyme fragment too short to be reliable)
    7. too large (one of the read-end mapped in a restriction-enzyme fragment too long to be reliable, e.g. may be centromeric)
    8. over-represented (not applied here)
    9. duplicated (both read ends mapped in the exact same way as an other read-end pair)
    10. random breaks (one ofthe read-end is mapped too far* from any possible restriction-enzyme cut site)
   
(*) *Too far*  and *too close* are defined automatically by TADbit inferring the distribution of sequenced fragment sizes from mapped dangling-ends reads.

The command following Nutlin1a-1h example would be:

```
tadbit filter -w HCT116_WT_NUT_1h -C 12 --noX  --apply 1 2 3 4 6 7 9 10 --clean --valid --compress_input
```

## 4. Normalizing

Interaction matrix with valid di-tags generatedin the previous step is now normalized. We used Vanilla ([Rao et al. 2014](https://doi.org/10.1016/j.cell.2014.11.021)) normalization that divides each cell of the matrix by the multiplication ofthe sum of interactions in its row and its column.

```
tadbit normalize -w HCT116_WT_NUT_1h -r 1000000 -C 24 --normalization Vanilla --renz HindIII --noX -F 1 2 3 4 6 7 9 10
```

## 5. Compartments calling

Compartments were called using individual chromosomic normalized matrices. These matrices were further normalized by distance, and a Pearson correlation matrix was then computed. First and second eigen vector of these matrices were kept. A  compartments were associated a given eigenvector sign (positive or negative) depending on its GC content.

The choice of which eigenvector to use (first or second) was decided by manual inspection of each result chromosome by chromosome (see Notebooks about `Compartments-parameter exploration`).

```
for chrom in {{1..22},X}
do
  tadbit segment -w HCT116_WT_NUT_1h --cpu 8 --reso 100000 --only_compartments --noX --fasta hg19.fa --reso 100000 -c $chrom --savecorr --n_evs 3 --ev_index 1 --format png --smoothing_window 3
done
```

## 6. TAD border calling

TADs were called at 50kb using TADbit internal methodology.

```
for chrom in {{1..22},X}
do
  tadbit segment -w HCT116_WT_NUT_1h --cpu 16 --only_tads --noX --reso 50000 -c $chrom
done
```

TAD borders generated are then aligned and insulation score along the genome is computed. See the `Align TAD-borders and compute insulation scores` notebooks.
