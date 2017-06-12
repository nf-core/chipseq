# ![NGI-ChIPseq](https://raw.githubusercontent.com/SciLifeLab/NGI-ChIPseq/master/docs/images/NGI-ChIPseq_logo.png)

[![Build Status](https://travis-ci.org/SciLifeLab/NGI-ChIPseq.svg?branch=master)](https://travis-ci.org/SciLifeLab/NGI-ChIPseq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg
)](https://www.nextflow.io/)

## Introduction
NGI-ChIPseq is a bioinformatics best-practice analysis pipeline used for ChIP-seq (chromatin immunoprecipitation sequencing) data analysis at the [National Genomics Infastructure](https://ngisweden.scilifelab.se/) at [SciLifeLab Stockholm](https://www.scilifelab.se/platforms/ngi/), Sweden.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

This pipeline is primarily used with a SLURM cluster on the Swedish [UPPMAX systems](https://www.uppmax.uu.se). However, the pipeline should be able to run on any system that Nextflow supports.

## Installation
### NextFlow installation
See https://github.com/SciLifeLab/NGI-NextflowDocs for instructions on how to install and configure
Nextflow.

### Pipeline installation
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub when run if
`SciLifeLab/NGI-ChIPseq` is specified as the pipeline name:
```bash
nextflow run SciLifeLab/NGI-ChIPseq
```

If you prefer, you can download the files yourself from GitHub and run them directly:
```bash
git clone https://github.com/SciLifeLab/NGI-ChIPseq.git
nextflow run NGI-ChIPseq/main.nf
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```
nextflow run SciLifeLab/NGI-ChIPseq --reads '*_R{1,2}.fastq.gz' --macsconfig 'macssetup.config'
```

Note that the pipeline will create files in your working directory:
```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### `--reads`
Location of the input FastQ files.

Default: `data/*{1,2}*.fastq.gz`

Single-end reads:
```
 --reads 'path/to/data/*.fastq.gz'
```

Paired-end reads:
```
 --reads 'path/to/data/sample_*_{1,2}.fastq.gz'
```

**NB: Must be enclosed in quotes!**

Note that the `{1,2}` parentheses are required to specify paired end data. Running `--reads '*.fastq'` will treat
all files as single end. The file path should be in quotation marks to prevent shell glob expansion.

If left unspecified, the pipeline will assume that the data is in a directory called `data` in the working directory.


### `--macsconfig`
The setup file for peak calling using MACS.

Default: `data/macsconfig`

Format:
```
ChIPSampleID1,CtrlSampleID1,AnalysisID1
ChIPSampleID2,CtrlSampleID2,AnalysisID2
ChIPSampleID3,,AnalysisID3
```

1. Column 1: ChIP sample name
  * Typically the file basename, shared between both reads
  * _eg._ for `chip_sample_1_R1.fastq.gz` and `chip_sample_1_R2.fastq.gz`, enter `chip_sample_1`
2. Column 2: Control sample name
  * Typically the input file basename, shared between both reads
  * _eg._ for `chip_input_R1.fastq.gz` and `chip_input_R2.fastq.gz`, enter `chip_input`
3. Column 3: Analysis ID
  * The analysis ID. Used for the output directory name.
  * Should be unique for each sample / line in the file.

_NB:_ For single-sample peaking calling without a control sample, skip the field of `CtrlSampleID`.

### `--genome`
The reference genome to use of the analysis, needs to be one of the genome specified in the config file.
Default: `None` (default can be set using a Nextflow configuration file, see below).

*Currently only the human and mouse genomes are fully supported by this pipeline. For other genomes `MACS` and `NGSplot` will not run!*
* Human
  * `--genome GRCh37`
* Mouse
  * `--genome GRCm38`

See [`conf/uppmax.config`](conf/uppmax.config) for a list of the reference genomes and their keys.

### `--extendReadsLen`
Amount of base pairs to extend the reads for the deepTools analysis.

Default: `100`

### `--bwa_index`
Instead of specifying genome and using the UPPMAX reference genome indices, you can specify the full
path to your reference genome when you run the pipeline:
```
--bwa_index [path to BWA index]
```

### `--outdir`
The output directory where the results will be saved.

Default: `./results`

### `--saveReference`
Supply this parameter to save any generated reference genome files to your results folder.
These can then be used for future pipeline runs, reducing processing times.

### `--saveTrimmed`
By default, trimmed FastQ files will not be saved to the results directory. Specify this
flag (or set to true in your config file) to copy these files when complete.

### `--saveAlignedIntermediates`
As above, by default intermediate BAM files will not be saved. The final BAM files created
after the Picard MarkDuplicates step are always saved. Set to true to also copy out BAM
files from BWA and sorting steps.

### `-c`
Specify the path to a specific config file (this is a core NextFlow command). Useful if using different UPPMAX
projects or different sets of reference genomes. **NB:** one hyphen only (core Nextflow parameter).


## Configuration
By default, the pipeline is configured to run on the Swedish UPPMAX cluster (milou / irma).

You will need to specify your UPPMAX project ID when running a pipeline. To do this, use
the command line flag `--project <project_ID>`.

To avoid having to specify this every time you run Nextflow, you can add it to your
personal Nextflow config file instead. Add this line to `~/.nextflow/config`:

```groovy
params.project = 'project_ID'
```

The pipeline will exit with an error message if you try to run it pipeline with the default
UPPMAX config profile and don't set project.

The same user config setup can be used to set a default genome, or any other command line parameters.
For example:

```groovy
params.genome = 'GRCh37'
params.extendReadsLen = 200
```


### Running on other clusters
It is entirely possible to run this pipeline on other clusters, though you will need to set up
your own config file so that the script knows where to find your reference files and how your
cluster works.

Copy the contents of [`conf/uppmax.config`](conf/uppmax.config) to your own config file somewhere
and then reference it with `-c` when running the pipeline.

If you think that there are other people using the pipeline who would benefit from your configuration
(eg. other common cluster setups), please let us know. It should be easy to create a new config file
in `conf` and reference this as a named profile in [`nextflow.config`](nextflow.config). Then these
configuration options can be used by specifying `-profile <name>` when running the pipeline.

## Credits
These scripts were written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/)
at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.
Written by Chuan Wang ([@chuan-wang](https://github.com/chuan-wang)) and Phil Ewels ([@ewels](https://github.com/ewels)).

<p align="center"><a href="stand_alone/http://www.scilifelab.se/" target="_blank"><img src="https://raw.githubusercontent.com/SciLifeLab/NGI-ChIPseq/master/docs/images/SciLifeLab_logo.png" title="SciLifeLab"></a></p>
