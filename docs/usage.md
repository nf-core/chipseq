# nf-core/chipseq: Usage

## Introduction

You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 6 columns, and a header row as shown in the examples below.

```bash
--input '[path to design file]'
```

### Multiple replicates

The `group` identifier should be identical when you have multiple replicates from the same experimental group, just increment the `replicate` identifier appropriately. The first replicate value for any given experimental group must be 1.

The `antibody` column is required to separate the downstream consensus peak merging and differential analysis for different antibodies. Its not advisable to generate a consensus peak set across different antibodies especially if their binding patterns are inherently different e.g. narrow transcription factors and broad histone marks.

The `control` column should be the `group` identifier for the controls for any given IP. The pipeline will automatically pair the inputs based on replicate identifier (i.e. where you have an equal number of replicates for your IP's and controls), alternatively, the first control sample in that group will be selected.

In the single-end design below there are triplicate samples for the `WT_BCATENIN_IP` group along with triplicate samples for their corresponding `WT_INPUT` samples.

```bash
group,replicate,fastq_1,fastq_2,antibody,control
WT_BCATENIN_IP,1,BLA203A1_S27_L006_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_BCATENIN_IP,2,BLA203A25_S16_L002_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_BCATENIN_IP,3,BLA203A49_S40_L001_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_INPUT,1,BLA203A6_S32_L006_R1_001.fastq.gz,,,
WT_INPUT,2,BLA203A30_S21_L002_R1_001.fastq.gz,,,
WT_INPUT,3,BLA203A31_S21_L003_R1_001.fastq.gz,,,
```

### Multiple runs of the same library

Both the `group` and `replicate` identifiers should be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will perform the alignments in parallel, and subsequently merge them before further analysis. Below is an example where the second replicate of the `WT_BCATENIN_IP` and `WT_INPUT` groups has been re-sequenced multiple times:

```bash
group,replicate,fastq_1,fastq_2,antibody,control
WT_BCATENIN_IP,1,BLA203A1_S27_L006_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_BCATENIN_IP,2,BLA203A25_S16_L001_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_BCATENIN_IP,2,BLA203A25_S16_L002_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_BCATENIN_IP,2,BLA203A25_S16_L003_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_BCATENIN_IP,3,BLA203A49_S40_L001_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_INPUT,1,BLA203A6_S32_L006_R1_001.fastq.gz,,,
WT_INPUT,2,BLA203A30_S21_L001_R1_001.fastq.gz,,,
WT_INPUT,2,BLA203A30_S21_L002_R1_001.fastq.gz,,,
WT_INPUT,3,BLA203A31_S21_L003_R1_001.fastq.gz,,,
```

### Full design

A final design file may look something like the one below. This is for two antibodies and associated controls in triplicate, where the second replicate of the `WT_BCATENIN_IP` and `NAIVE_BCATENIN_IP` group has been sequenced twice:

```bash
group,replicate,fastq_1,fastq_2,antibody,control
WT_BCATENIN_IP,1,BLA203A1_S27_L006_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_BCATENIN_IP,2,BLA203A25_S16_L001_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_BCATENIN_IP,2,BLA203A25_S16_L002_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_BCATENIN_IP,3,BLA203A49_S40_L001_R1_001.fastq.gz,,BCATENIN,WT_INPUT
NAIVE_BCATENIN_IP,1,BLA203A7_S60_L001_R1_001.fastq.gz,,BCATENIN,NAIVE_INPUT
NAIVE_BCATENIN_IP,2,BLA203A43_S34_L001_R1_001.fastq.gz,,BCATENIN,NAIVE_INPUT
NAIVE_BCATENIN_IP,2,BLA203A43_S34_L002_R1_001.fastq.gz,,BCATENIN,NAIVE_INPUT
NAIVE_BCATENIN_IP,3,BLA203A64_S55_L001_R1_001.fastq.gz,,BCATENIN,NAIVE_INPUT
WT_TCF4_IP,1,BLA203A3_S29_L006_R1_001.fastq.gz,,TCF4,WT_INPUT
WT_TCF4_IP,2,BLA203A27_S18_L001_R1_001.fastq.gz,,TCF4,WT_INPUT
WT_TCF4_IP,3,BLA203A51_S42_L001_R1_001.fastq.gz,,TCF4,WT_INPUT
NAIVE_TCF4_IP,1,BLA203A9_S62_L001_R1_001.fastq.gz,,TCF4,NAIVE_INPUT
NAIVE_TCF4_IP,2,BLA203A45_S36_L001_R1_001.fastq.gz,,TCF4,NAIVE_INPUT
NAIVE_TCF4_IP,3,BLA203A66_S57_L001_R1_001.fastq.gz,,TCF4,NAIVE_INPUT
WT_INPUT,1,BLA203A6_S32_L006_R1_001.fastq.gz,,,
WT_INPUT,2,BLA203A30_S21_L001_R1_001.fastq.gz,,,
WT_INPUT,3,BLA203A31_S21_L003_R1_001.fastq.gz,,,
NAIVE_INPUT,1,BLA203A12_S3_L001_R1_001.fastq.gz,,,
NAIVE_INPUT,2,BLA203A48_S39_L001_R1_001.fastq.gz,,,
NAIVE_INPUT,3,BLA203A49_S1_L006_R1_001.fastq.gz,,,
```

| Column      | Description                                                                                                                                      |
|-------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| `group`     | Group/condition identifier for sample. This will be identical for re-sequenced libraries and replicate samples from the same experimental group. |
| `replicate` | Integer representing replicate number. This will be identical for re-sequenced libraries. Must start from `1..<number of replicates>`.           |
| `fastq_1`   | Full path to FastQ file for read 1. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".                                        |
| `fastq_2`   | Full path to FastQ file for read 2. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".                                        |
| `antibody`  | Antibody name. This is required to segregate downstream analysis for different antibodies. Required when `control` is specified.                 |
| `control`   | Group identifier for control sample. The pipeline will automatically select the control sample with the same replicate identifier as the IP.     |

Example design files have been provided with the pipeline for [paired-end](../assets/design_pe.csv) and [single-end](../assets/design_se.csv) data.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/chipseq --input design.csv --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/chipseq
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/chipseq releases page](https://github.com/nf-core/chipseq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
    * A generic configuration profile to be used with [Docker](https://docker.com/)
    * Pulls software from Docker Hub: [`nfcore/chipseq`](https://hub.docker.com/r/nfcore/chipseq/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
    * Pulls software from Docker Hub: [`nfcore/chipseq`](https://hub.docker.com/r/nfcore/chipseq/)
* `conda`
    * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
    * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
    * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters

### `--input`

You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 6 columns, and a header row as shown in the examples below.

```bash
--input '[path to design file]'
```

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core NextFlow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `star` 32GB of memory, you could use the following config:

```nextflow
process {
  withName: star {
    memory = 32.GB
  }
}
```

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

#### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
