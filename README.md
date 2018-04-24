# ![nf-core/ChIPseq](docs/images/chipseq_logo.png)

[![Build Status](https://travis-ci.org/nf-core/ChIPseq.svg?branch=master)](https://travis-ci.org/nf-core/ChIPseq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.6-brightgreen.svg)](https://www.nextflow.io/)
[![Gitter](https://img.shields.io/badge/gitter-%20join%20chat%20%E2%86%92-4fb99a.svg)](https://gitter.im/nf-core/Lobby)

[![Docker Container available](https://img.shields.io/docker/automated/nfcore/chipseq.svg)](https://hub.docker.com/r/nfcore/chipseq/)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

## Introduction
**nf-core/ChIPseq** is a bioinformatics best-practice analysis pipeline used for chromatin immunoprecipitation (ChIP-seq) data analysis.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

### Pipeline Steps

* Make BWA reference genome index (optional)
* FastQC
* TrimGalore!
* Align with BWA
* Samtools sort, index, stats & convert to BED
* Samtools idxstats
* Picard MarkDuplicates
* Count read statistics
* PhantomPeakQualTools
* SPP: Calculate NSCRSC and cross correlation
* DeepTools bamPEFragmentSize, plotFingerprint, bamCoverage, multiBamSummary, plotCorrelation, plotPCA
* NGSplot
* MACS2 peak calling
* Saturation analysis
* ChIP peak annotation
* MultiQC


### Documentation
The nf-core/ChIPseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation and configuration](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Output and how to interpret the results](docs/output.md)


## Credits
These scripts were written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/)
at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.
Written by Chuan Wang ([@chuan-wang](https://github.com/chuan-wang)) and Phil Ewels ([@ewels](https://github.com/ewels)).
