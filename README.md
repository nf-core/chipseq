# ![nf-core/chipseq](docs/images/nfcore-chipseq_logo.png)

[![Build Status](https://travis-ci.org/nf-core/chipseq.svg?branch=master)](https://travis-ci.org/nf-core/chipseq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/chipseq.svg)](https://hub.docker.com/r/nfcore/chipseq/)

## Introduction
**nf-core/chipseq** is a bioinformatics best-practice analysis pipeline used for chromatin immunoprecipitation (ChIP-seq) data analysis.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

### Pipeline Steps

* Make BWA reference genome index (optional)
* FastQC raw sequencing reads quality control
* TrimGalore! adapter trimming
* BWA alignment against reference genome
* Samtools post-alignment processing and format conversion
* Picard duplicate read identification
* Statistics about read counts
* PhantomPeakQualTools normalized strand cross-correlation (NSC) and relative strand cross-correlation (RSC)
* deepTools fingerprint, coverage bigwig, correlation plots of reads over genome-wide bins, and distribution of reads around gene bodies
* MACS2 peak calling
* Saturation analysis (optional)
* Post peak calling processing: blacklist filtering and annotation
* MultiQC

### Documentation
The nf-core/chipseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
    * [Reference genomes](docs/configuration/reference_genomes.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

## Credits
These scripts were written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/)
at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.
Written by Chuan Wang ([@chuan-wang](https://github.com/chuan-wang)) and Phil Ewels ([@ewels](https://github.com/ewels)).
