# nf-core/ChIPseq Installation

To start using the nf-core/ChIPseq pipeline, there are three steps described below:

1. [Install Nextflow](#1-install-nextflow)
2. [Install the pipeline](#2-install-the-pipeline)
3. Configure the pipeline
    * [Swedish UPPMAX System](#31-configuration-uppmax)
    * [Other Clusters](#32-configuration-other-clusters)
    * [Docker](#33-configuration-docker)
    * [Amazon AWS](#34-configuration-amazon-ec2)

## 1) Install NextFlow
Nextflow runs on most POSIX systems (Linux, Mac OSX etc). It can be installed by running the following commands:

```bash
# Make sure that Java v7+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your PATH:
mv nextflow ~/bin
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

## 2) Install the Pipeline
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub if `nf-core/ChIPseq` is specified as the pipeline name.

### Offline use

If you need to run the pipeline on a system with no internet connection, you will need to download the files yourself from GitHub and run them directly:

```bash
wget https://github.com/nf-core/ChIPseq/archive/master.zip
unzip master.zip -d /my-pipelines/
cd /my_data/
nextflow run /my-pipelines/nf-core/ChIPseq-master
```

## 3.1) Configuration: UPPMAX
The pipeline comes bundled with configurations to use the [Swedish UPPMAX](https://www.uppmax.uu.se/) clusters (tested on `milou`, `rackham`, `bianca` and `irma`). As such, you shouldn't need to add any custom configuration - everything _should_ work out of the box.

To use the pipeline on UPPMAX, you **must** specificy `-profile uppmax` when running the pipeline (as of Nov 2017).

Note that you will need to specify your UPPMAX project ID when running a pipeline. To do this, use the command line flag `--project <project_ID>`. The pipeline will exit with an error message if you try to run it pipeline with the UPPMAX config profile without a project.

**Optional Extra:** To avoid having to specify your project every time you run Nextflow, you can add it to your personal Nextflow config file instead. Add this line to `~/.nextflow/config`:

```nextflow
params.project = 'project_ID' // eg. b2017123
```

## 3.2) Configuration: Other clusters
It is entirely possible to run this pipeline on other clusters, though you will need to set up your own config file so that the script knows where to find your reference files and how your cluster works.

> If you think that there are other people using the pipeline who would benefit from your configuration (eg. other common cluster setups), please let us know. We can add a new configuration and profile which can used by specifying `-profile <name>` when running the pipeline.

If you are the only person to be running this pipeline, you can create your config file as `~/.nextflow/config` and it will be applied every time you run Nextflow. Alternatively, save the file anywhere and reference it when running the pipeline with `-c path/to/config`.

A basic configuration comes with the pipeline, which runs by default (the `standard` config profile with [`conf/base.config`](../conf/base.config)). This means that you only need to configure the specifics for your system and overwrite any defaults that you want to change.

### Cluster Environment
By default, Nextflow uses the `local` executor - in other words, all jobs are run in the login session. If you're using a simple server, this may be fine. If you're using a compute cluster, this is bad as all jobs will run on the head node.

To specify your cluster environment, add the following line to your config file:

```nextflow
process {
  executor = 'YOUR_SYSTEM_TYPE'
}
```

Many different cluster types are supported by Nextflow. For more information, please see the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html).

Note that you may need to specify cluster options, such as a project or queue. To do so, use the `clusterOptions` config option:

```nextflow
process {
  executor = 'SLURM'
  clusterOptions = '-A myproject'
}
```

### Reference Genomes
The nf-core/ChIPseq pipeline needs a reference genome for alignment and annotation. If not already available, start by downloading the relevant reference, for example from [illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html).

A reference genome path can be specified on the command line each time you run with `--bwa_index`. If no BWA index
is available, one can be generated using a FASTA file supplied with `--fasta`.
Alternatively, add the paths to the config under a relevant id and just specify this id with `--genome ID` when you run the pipeline _(this can also be set as a default in your config)_:

```nextflow
params {
  genomes {
    'YOUR-ID' {
      bwa  = '<path to the bwa index folder>'
      fasta = '<path to the fasta file>' // used if bwa index not given
    }
    'OTHER-GENOME' {
      // [..]
    }
  }
  // Optional - default genome. Ignored if --genome 'OTHER-GENOME' specified on command line
  genome = 'YOUR-ID'
}
```


### Software Requirements
To run the pipeline, several software packages are required. How you satisfy these requirements is essentially up to you and depends on your system.

#### Docker Image
Docker is a software tool that allows you to run your analysis inside a "container" - basically a miniature self-contained software environment. We've already made a docker image for you, so if you can run docker and nextflow then you don't need to worry about any other software dependencies.

The pipeline comes with a script to build a docker image. This runs automatically and creates a hosted docker image that you can find here: https://hub.docker.com/r/nf-core/ChIPseq/

If you run the pipeline with `-profile docker` or `-with-docker 'nf-core/ChIPseq'` then nextflow will download this docker image automatically and run using this.

Note that the docker images are tagged with version as well as the code, so this is a great way to ensure reproducibility. You can specify pipeline version when running with `-r`, for example `-r v1.4`. This uses pipeline code and docker image from this tagged version.

#### Singularity image
Many HPC environments are not able to run Docker due to problems with needing administrator privileges. [Singularity](http://singularity.lbl.gov/) is a tool designed to run on such HPC systems which is very similar to Docker. Even better, it can use create images directly from dockerhub. The UPPMAX configuration uses Singularity by default, meaning no problems with software dependencies and great reproducibility.

To use the singularity image with a different config, use `-with-singularity 'docker://nf-core/ChIPseq'` when running the pipeline.

If you intend to run the pipeline offline, nextflow will not be able to automatically download the singularity image for you. Instead, you'll have to do this yourself manually first, transfer the image file and then point to that.

First, pull the image file where you have an internet connection:

```bash
singularity pull --name nfcore-chipseq.img docker://nf-core/ChIPseq
```

Then transfer this file and run the pipeline with this path:

```bash
nextflow run /path/to/nf-core/ChIPseq -with-singularity /path/to/nfcore-chipseq.img
```

#### Environment Modules
If your cluster uses environment modules, you can use the pipeline with these. There is a bundled config file to use these on UPPMAX (as was done in earlier versions of this pipeline). To use this, run the pipeline with `-profile uppmax_modules`.

If running on another system, add lines to your custom config file as follows _(customise module names and versions as appropriate)_:

```nextflow
process {
  $fastqc.module = ['FastQC']
  $trim_galore.module = ['FastQC', 'TrimGalore']
  $bw.module = ['bwa', 'samtools/1.3']
  $samtools.module = ['samtools/1.3', 'BEDTools']
  $picard.module = ['picard/2.0.1', 'samtools/1.3', 'BEDTools']
  $phantompeakqualtools.module = ['R/3.2.3', 'phantompeakqualtools']
  $deepTools.module = ['deepTools']
  $ngsplot.module = ['samtools/1.3', 'R/3.2.3', 'ngsplot']
  $macs.module = ['MACS', 'samtools/1.3']
  $multiqc.module = ['MultiQC']
}
```

#### Manual Installation
If the software is not already available, you will need to install it.

If you are able to use [Docker](https://www.docker.com/), you can use the [sclifelab/nfcore-chipseq](https://hub.docker.com/r/nf-core/ChIPseq/) image which comes with all requirements. This is pulled by Nextflow automatically if you use `-profile docker` (see below for [further instructions](#33-configuration-docker)).

We recommend using [Bioconda](https://bioconda.github.io/) to install the required software as the process is quite easy in our experience.

## 3.3) Configuration: Docker
Docker is a great way to run nf-core/ChIPseq, as it manages all software installations and allows the pipeline to be run in an identical software environment across a range of systems.

Nextflow has [excellent integration](https://www.nextflow.io/docs/latest/docker.html) with Docker, and beyond installing the two tools, not much else is required.

First, install docker on your system : [Docker Installation Instructions](https://docs.docker.com/engine/installation/)

Then, simply run the analysis pipeline:
```bash
nextflow run nf-core/ChIPseq -profile docker # rest of normal launch command
```

Nextflow will recognise `nf-core/ChIPseq` and download the pipeline from GitHub. The `-profile docker` configuration lists the [sclifelab/nfcore-chipseq](https://hub.docker.com/r/nf-core/ChIPseq/) image that we have created and is hosted at dockerhub, and this is downloaded.

A reference genome is still required by the pipeline. See the above [Reference Genomes](#reference-genomes) documentation for instructions on how to configure Nextflow with preset paths to make this easier.

A test suite for docker comes with the pipeline, and can be run by moving to the [`tests` directory](https://github.com/ewels/nf-core/ChIPseq/tree/master/tests) and running `./docker_test.sh`. This will download a small lambda genome and some data, and attempt to run the pipeline through docker on that small dataset. This is automatically run using [Travis](https://travis-ci.org/nf-core/ChIPseq/) whenever changes are made to the pipeline.
