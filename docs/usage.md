# nf-core/chipseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/chipseq/usage](https://nf-co.re/chipseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 5 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple replicates

The `sample` identifier should be identical when you have multiple replicates from the same experimental group, just increment the `replicate` identifier appropriately. The first replicate value for any given experimental group must be 1.

The `antibody` column is required to separate the downstream consensus peak merging for different antibodies. It is not advisable to generate a consensus peak set across different antibodies especially if their binding patterns are inherently different e.g. narrow transcription factors and broad histone marks.

The `control` column should be the `sample` identifier for the controls for any given IP. This column together with the `control_replicate` column will set the corresponding control for each of the samples in the table.

```console
group,fastq_1,fastq_2,replicate,antibody,control,control_replicate
WT_BCATENIN_IP,BLA203A1_S27_L006_R1_001.fastq.gz,,1,BCATENIN,WT_INPUT,1
WT_BCATENIN_IP,BLA203A25_S16_L002_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
WT_BCATENIN_IP,BLA203A49_S40_L001_R1_001.fastq.gz,,3,BCATENIN,WT_INPUT,3
WT_INPUT,BLA203A6_S32_L006_R1_001.fastq.gz,,1,,,
WT_INPUT,BLA203A30_S21_L002_R1_001.fastq.gz,,2,,,
WT_INPUT,BLA203A31_S21_L003_R1_001.fastq.gz,,3,,,
```

### Multiple runs of the same library

Both the `sample` and `replicate` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will perform the alignments in parallel, and subsequently merge them before further analysis. Below is an example where the samples called `WT_BCATENIN_IP` and `WT_INPUT` have been re-sequenced multiple times:

```console
sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate
WT_BCATENIN_IP,BLA203A1_S27_L006_R1_001.fastq.gz,,1,BCATENIN,WT_INPUT,1
WT_BCATENIN_IP,BLA203A25_S16_L001_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
WT_BCATENIN_IP,BLA203A25_S16_L002_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
WT_BCATENIN_IP,BLA203A25_S16_L003_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
WT_BCATENIN_IP,BLA203A49_S40_L001_R1_001.fastq.gz,,3,BCATENIN,WT_INPUT,3
WT_INPUT,BLA203A6_S32_L006_R1_001.fastq.gz,,1,,,
WT_INPUT,BLA203A30_S21_L001_R1_001.fastq.gz,,2,,,
WT_INPUT,BLA203A30_S21_L002_R1_001.fastq.gz,,2,,,
WT_INPUT,BLA203A31_S21_L003_R1_001.fastq.gz,,3,,,
```

### Full design

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 7 columns to match those defined in the table below.

A final design file may look something like the one below. This is for two antibodies and associated controls, where the second replicate of the `WT_BCATENIN_IP` and `NAIVE_BCATENIN_IP` samples have been sequenced twice:

```console
sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate
WT_BCATENIN_IP,BLA203A1_S27_L006_R1_001.fastq.gz,,1,BCATENIN,WT_INPUT,1
WT_BCATENIN_IP,BLA203A25_S16_L001_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
WT_BCATENIN_IP,BLA203A25_S16_L002_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
WT_BCATENIN_IP,BLA203A49_S40_L001_R1_001.fastq.gz,,3,BCATENIN,WT_INPUT,3
NAIVE_BCATENIN_IP,BLA203A7_S60_L001_R1_001.fastq.gz,,1,BCATENIN,NAIVE_INPUT,1
NAIVE_BCATENIN_IP,BLA203A43_S34_L001_R1_001.fastq.gz,,2,BCATENIN,NAIVE_INPUT,2
NAIVE_BCATENIN_IP,BLA203A43_S34_L002_R1_001.fastq.gz,,2,BCATENIN,NAIVE_INPUT,2
NAIVE_BCATENIN_IP,BLA203A64_S55_L001_R1_001.fastq.gz,,3,BCATENIN,NAIVE_INPUT,3
WT_TCF4_IP,BLA203A3_S29_L006_R1_001.fastq.gz,,1,TCF4,WT_INPUT,1
WT_TCF4_IP,BLA203A27_S18_L001_R1_001.fastq.gz,,2,TCF4,WT_INPUT,2
WT_TCF4_IP,BLA203A51_S42_L001_R1_001.fastq.gz,,3,TCF4,WT_INPUT,3
NAIVE_TCF4_IP,BLA203A9_S62_L001_R1_001.fastq.gz,,1,TCF4,NAIVE_INPUT,1
NAIVE_TCF4_IP,BLA203A45_S36_L001_R1_001.fastq.gz,,2,TCF4,NAIVE_INPUT,2
NAIVE_TCF4_IP,BLA203A66_S57_L001_R1_001.fastq.gz,,3,TCF4,NAIVE_INPUT,3
WT_INPUT,BLA203A6_S32_L006_R1_001.fastq.gz,,1,,,
WT_INPUT,BLA203A30_S21_L001_R1_001.fastq.gz,,2,,,
WT_INPUT,BLA203A31_S21_L003_R1_001.fastq.gz,,3,,,
NAIVE_INPUT,BLA203A12_S3_L001_R1_001.fastq.gz,,1,,,
NAIVE_INPUT,BLA203A48_S39_L001_R1_001.fastq.gz,,2,,,
NAIVE_INPUT,BLA203A49_S1_L006_R1_001.fastq.gz,,3,,,
```

| Column              | Description                                                                                                                                                                            |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`            | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`           | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`           | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `replicate`         | Integer representing replicate number. This will be identical for re-sequenced libraries. Must start from `1..<number of replicates>`.                                                 |
| `antibody`          | Antibody name. This is required to segregate downstream analysis for different antibodies. Required when `control` is specified.                                                       |
| `control`           | Sample name for control sample.                                                                                                                                                        |
| `control_replicate` | Integer representing replicate number for control sample.                                                                                                                              |

Example design files have bee_n provided with the pipeline for [paired-end](../assets/samplesheet_pe.csv) and [single-end](../assets/samplesheet_se.csv) data.

> **NB:** The `group` and `replicate` columns were replaced with a single `sample` column as of v2.0 of the pipeline. The `sample` column is essentially a concatenation of the `group` and `replicate` columns. If all values of `sample` have the same number of underscores, fields defined by these underscore-separated names may be used in the PCA plots produced by the pipeline, to regain the ability to represent different groupings.

## Reference genome files

The minimum reference genome requirements are a FASTA and GTF file, all other files required to run the pipeline can be generated from these files. However, it is more storage and compute friendly if you are able to re-use reference genome files as efficiently as possible. It is recommended to use the `--save_reference` parameter if you are using the pipeline to build new indices (e.g. those unavailable on [AWS iGenomes](https://nf-co.re/usage/reference_genomes)) so that you can save them somewhere locally. The index building step can be quite a time-consuming process and it permits their reuse for future runs of the pipeline to save disk space. You can then either provide the appropriate reference genome files on the command-line via the appropriate parameters (e.g. `--bwa_index '/path/to/bwa/index/'`) or via a custom config file.

- If `--genome` is provided then the FASTA and GTF files (and existing indices) will be automatically obtained from AWS-iGenomes unless these have already been downloaded locally in the path specified by `--igenomes_base`.
- If `--gene_bed` is not provided then it will be generated from the GTF file.

> **NB:** Compressed reference files are also supported by the pipeline i.e. standard files with the `.gz` extension and indices folders with the `tar.gz` extension.

## Blacklist bed files

The blacklist bed files where obtained using the commands below:

```console
cd ..
mkdir -p v1.0
cd v1.0
wget -L https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz && gunzip ENCFF001TDO.bed.gz && mv ENCFF001TDO.bed hg19-blacklist.v1.bed

mkdir -p assets/blacklists/v2.0/
cd assets/blacklists/v2.0/
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce10-blacklist.v2.bed.gz && gunzip ce10-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce11-blacklist.v2.bed.gz && gunzip ce11-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm3-blacklist.v2.bed.gz && gunzip dm3-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm6-blacklist.v2.bed.gz && gunzip dm6-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz && gunzip hg19-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz && gunzip hg38-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz && gunzip mm10-blacklist.v2.bed.gz

cd ..
mkdir -p v3.0
cd v3.0
wget -L https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz && gunzip ENCFF356LFX.bed.gz && mv ENCFF356LFX.bed hg38-blacklist.v3.bed
```

> **NB:** A detailed description of the different versions of the files can be found [here](https://sites.google.com/site/anshulkundaje/projects/blacklists). Also, to to see which blacklist bed files are assigned by default to the respective reference genome check the [igenomes.config](https://github.com/nf-core/chipseq/blob/master/conf/igenomes.config).

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/chipseq --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> ⚠️ Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/chipseq -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/chipseq
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/chipseq releases page](https://github.com/nf-core/chipseq/releases) and find the latest pipeline version - numeric only (eg. `2.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 2.0.0`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> 💡 If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
