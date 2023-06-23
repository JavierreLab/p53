# ChIP-seq processing

Here, you can find a detailed description of how ChIP-seq data was processed in this project.

### Dependencies

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


### Summary of the workflow

1. **Trimming and Quality**: Trim Galore and FastQC
2. **Alignment**: Bowtie2
3. **Filtering**: samtools and sambamba
4. **Coverage**: deepTools
5. **Peak Calling**: Macs2
6. **Generating Peakmatrix**: csaw and GenomicRanges
7. **Background Correction**: csaw and GenomicRanges
8. **Differential Analysis**: DESeq2


A detailed description of how steps 1-5 were performed can be found [here](https://github.com/JavierreLab/liCHiC/tree/main/3.ChIPseq%20Processing), the github page associated to Tomás-Daza, L *et al.* Low input capture Hi-C (liCHi-C) identifies promoter-enhancer interactions at high-resolution. *Nature Communications* **14**, 268 (2023). [https://doi.org/10.1038/s41467-023-35911-8](https://doi.org/10.1038/s41467-023-35911-8)

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

