# Hi-C processing

*This README will be finished soon* 

Here, you can find a detailed description of how RNA-seq data was processed in this project.

## Dependencies

* [TADbit](https://github.com/fransua/tadbit/tree/p53_javierre)

## Summary of the workflow

1. **Alignment**: hisat2 
2. **Parsing**: pull sequencing lanes and replicates together
3. **Filtering**: remove non-canonical di-tags
4. **Normalizing**: remove Hi-C specific biases
5. **Compartments calling**: call for A/B compartments
6. **TAD calling**: call for Topologically Associating Domains (TADs)


## 1. Alignment
