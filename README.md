# ![NGI-ChIPseq](https://raw.githubusercontent.com/SciLifeLab/NGI-ChIPseq/master/docs/images/NGI-ChIPseq_logo.png)

[![Build Status](https://travis-ci.org/SciLifeLab/NGI-ChIPseq.svg?branch=master)](https://travis-ci.org/SciLifeLab/NGI-ChIPseq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg
)](https://www.nextflow.io/)

## Introduction
NGI-ChIPseq is a bioinformatics best-practice analysis pipeline used for ChIP-seq (chromatin immunoprecipitation sequencing) data analysis at the [National Genomics Infastructure](https://ngisweden.scilifelab.se/) at [SciLifeLab Stockholm](https://www.scilifelab.se/platforms/ngi/), Sweden.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

This pipeline is primarily used with a SLURM cluster on the Swedish [UPPMAX systems](https://www.uppmax.uu.se). However, the pipeline should be able to run on any system that Nextflow supports. We have done some limited testing using Docker and AWS, and the pipeline comes with some configuration for these systems. See the [installation docs](docs/installation.md) for more information.

### Documentation
The NGI-ChIPseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation and configuration](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Output and how to interpret the results](docs/output.md)

If you're interested in running the pipeline in the cloud, please read the docs about using our pipeline with Amazon Web Services on the [NGI-ChIPseq pipeline](https://github.com/SciLifeLab/NGI-ChIPseq/blob/master/docs/amazon_web_services.md) (the instructions should work with this pipeline as well).

## Credits
These scripts were written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/)
at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.
Written by Chuan Wang ([@chuan-wang](https://github.com/chuan-wang)) and Phil Ewels ([@ewels](https://github.com/ewels)).

---

[![SciLifeLab](https://raw.githubusercontent.com/SciLifeLab/NGI-ChIPseq/master/docs/images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](https://raw.githubusercontent.com/SciLifeLab/NGI-ChIPseq/master/docs/images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---
