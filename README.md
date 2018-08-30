# ![nf-core/chipseq](docs/images/chipseq_logo.png)

[![Build Status](https://travis-ci.org/nf-core/chipseq.svg?branch=master)](https://travis-ci.org/nf-core/chipseq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.31.1-brightgreen.svg)](https://www.nextflow.io/)
[![Gitter](https://img.shields.io/badge/gitter-%20join%20chat%20%E2%86%92-4fb99a.svg)](https://gitter.im/nf-core/Lobby)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker Container available](https://img.shields.io/docker/automated/nfcore/chipseq.svg)](https://hub.docker.com/r/nfcore/chipseq/)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/1123)


## Introduction
**nf-core/chipseq** is a bioinformatics best-practice analysis pipeline used for chromatin immunoprecipitation (ChIP-seq) and assay for transposase accessible chromatin (ATAC-seq) data analysis.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

### Pipeline Steps

* Make BWA reference genome index (optional)
* Build BED reference based on GTF (optional)
* Build genome size table for bedToBam conversion (optional)
* FastQC for initial quality control of sequence reads
* TrimGalore! for adapter trimming
* BWA for alignment
* Samtools for post-alignment processing with and alignment statistics
* Picard MarkDuplicates for duplicate removal
* Count read statistics
* Phantompeakqualtools for NSC, RSC and strand-shift cross correlation plot
* DeepTools for paired-end fragment size distribution, fingerprint, reads distribution profile, sample pair-wise correlation, and PCA plot.
* MACS2 for peak calling
* MACS2 for saturation analysis (optional)
* Bioconductor ChIPpeakAnno for peak annotation
* MultiQC


### Documentation
The nf-core/chipseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation and configuration](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Output and how to interpret the results](docs/output.md)


## Credits
These scripts were written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/)
at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.
Written by Chuan Wang ([@chuan-wang](https://github.com/chuan-wang)) and Phil Ewels ([@ewels](https://github.com/ewels)).
