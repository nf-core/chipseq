# NGI-ChIPseq Installation

To start using the NGI-ChIPseq pipeline, there are three steps described below:

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

See [nextflow.io](https://www.nextflow.io/) and [NGI-NextflowDocs](https://github.com/SciLifeLab/NGI-NextflowDocs) for further instructions on how to install and configure Nextflow.

## 2) Install the Pipeline
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub if `SciLifeLab/NGI-ChIPseq` is specified as the pipeline name.

If you prefer, you can download the files yourself from GitHub and run them directly:

```bash
git clone https://github.com/SciLifeLab/NGI-ChIPseq.git
nextflow run NGI-ChIPseq/main.nf
```

## 3.1) Configuration: UPPMAX
By default, the pipeline is configured to run on the [Swedish UPPMAX](https://www.uppmax.uu.se/) cluster (`milou` / `irma`). As such, you shouldn't need to add any custom configuration - everything _should_ work out of the box.

Note that you will need to specify your UPPMAX project ID when running a pipeline. To do this, use the command line flag `--project <project_ID>`. The pipeline will exit with an error message if you try to run it pipeline with the default UPPMAX config profile without a project.

**Optional Extra:** To avoid having to specify your project every time you run Nextflow, you can add it to your personal Nextflow config file instead. Add this line to `~/.nextflow/config`:

```groovy
params.project = 'project_ID' // eg. b2017123
```

## 3.2) Configuration: Other clusters
It is entirely possible to run this pipeline on other clusters, though you will need to set up your own config file so that the script knows where to find your reference files and how your cluster works.

> If you think that there are other people using the pipeline who would benefit from your configuration (eg. other common cluster setups), please let us know. We can add a new configuration and profile which can used by specifying `-profile <name>` when running the pipeline.

If you are the only person to be running this pipeline, you can create your config file as `~/.nextflow/config` and it will be applied every time you run Nextflow. Alternatively, save the file anywhere and reference it when running the pipeline with `-c path/to/config`.

An empty configuration comes with the pipeline, which should be applied by using the command line flag `-profile none`. This prevents the UPPMAX defaults (above) from being applied and means that you only need to configure the specifics for your system.

### Cluster Environment
By default, Nextflow uses the `local` executor - in other words, all jobs are run in the login session. If you're using a simple server, this may be fine. If you're using a compute cluster, this is bad as all jobs will run on the head node.

To specify your cluster environment, add the following line to your config file:

```groovy
process {
  executor = 'YOUR_SYSTEM_TYPE'
}
```

Many different cluster types are supported by Nextflow. For more information, please see the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html).

Note that you may need to specify cluster options, such as a project or queue. To do so, use the `clusterOptions` config option:

```groovy
process {
  executor = 'SLURM'
  clusterOptions = '-A myproject'
}
```

### Reference Genomes
The NGI-ChIPseq pipeline needs a reference genome for alignment and annotation. If not already available, start by downloading the relevant reference, for example from [illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html).

A reference genome path can be specified on the command line each time you run with `--bwa_index`. If no BWA index
is available, one can be generated using a FASTA file supplied with `--fasta`.
Alternatively, add the paths to the config under a relevant id and just specify this id with `--genome ID` when you run the pipeline _(this can also be set as a default in your config)_:

```groovy
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

#### Environment Modules
If your cluster uses environment modules, the software may already be available. If so, just add lines to your custom config file as follows _(customise module names and versions as appropriate)_:

```groovy
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

If you are able to use [Docker](https://www.docker.com/), you can use the [sclifelab/ngi-chipseq](https://hub.docker.com/r/scilifelab/ngi-chipseq/) image which comes with all requirements. This is pulled by Nextflow automatically if you use `-profile docker` (see below for [further instructions](#33-configuration-docker)).

We recommend using [Bioconda](https://bioconda.github.io/) to install the required software as the process is quite easy in our experience.

## 3.3) Configuration: Docker
Docker is a great way to run NGI-ChIPseq, as it manages all software installations and allows the pipeline to be run in an identical software environment across a range of systems.

Nextflow has [excellent integration](https://www.nextflow.io/docs/latest/docker.html) with Docker, and beyond installing the two tools, not much else is required.

First, install docker on your system : [Docker Installation Instructions](https://docs.docker.com/engine/installation/)

Then, simply run the analysis pipeline:
```bash
nextflow run SciLifeLab/NGI-ChIPseq -profile docker # rest of normal launch command
```

Nextflow will recognise `SciLifeLab/NGI-ChIPseq` and download the pipeline from GitHub. The `-profile docker` configuration lists the [sclifelab/ngi-chipseq](https://hub.docker.com/r/scilifelab/ngi-chipseq/) image that we have created and is hosted at dockerhub, and this is downloaded.

A reference genome is still required by the pipeline. See the above [Reference Genomes](#reference-genomes) documentation for instructions on how to configure Nextflow with preset paths to make this easier.

A test suite for docker comes with the pipeline, and can be run by moving to the [`tests` directory](https://github.com/ewels/NGI-ChIPseq/tree/master/tests) and running `./docker_test.sh`. This will download a small lambda genome and some data, and attempt to run the pipeline through docker on that small dataset. This is automatically run using [Travis](https://travis-ci.org/SciLifeLab/NGI-ChIPseq/) whenever changes are made to the pipeline.

## 3.4) Configuration: Amazon EC2
There are multiple ways of running this pipeline over Amazon's EC2 service. Please see the [NGI-ChIPseq pipeline docs](https://github.com/SciLifeLab/NGI-ChIPseq/blob/master/docs/amazon_web_services.md) for more information.

---

[![SciLifeLab](images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---