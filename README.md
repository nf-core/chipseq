# NGI-ChIPseq
Nextflow ChIP-seq data analysis pipeline, National Genomics Infrastructure, Science for Life Laboratory in Stockholm, Sweden

Authors: Chuan Wang (@chuan-wang)

## Installation

##### NextFlow installation

To use this pipeline, you need to have a working version of NextFlow installed. You can find more information about this pipeline tool at [nextflow.io](http://www.nextflow.io/). The typical installation of NextFlow looks like this:

  `curl -fsSL get.nextflow.io | bash` 

  `mv ./nextflow ~/bin`

Note that if you're running on the Swedish UPPMAX cluster (Milou) you can load NextFlow as an environment module:

  `module load nextflow`

  `NextFlow configuration`
  
  
Next, you need to set up a config file so that NextFlow knows how to run and where to find reference indexes. You can find an example configuration file for UPPMAX (milou) with this repository: example_uppmax_config.

Copy this file to ~/.nextflow/config and edit the line '-A b2013064' to contain your own UPPMAX project identifier instead.

It is entirely possible to run this pipeline on other clusters - just note that you may need to customise the process environment (eg. if you're using a cluster system other than SLURM) and the paths to reference files.

##### Pipeline installation

This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub when run if SciLifeLab/NGI-ChIPseq is specified as the pipeline name.

If you prefer, you can download the files yourself from GitHub and run them directly:

`git clone https://github.com/SciLifeLab/NGI-ChIPseq.git`

`nextflow NGI-ChIPseq/main.nf`

## Running the pipeline

The typical command for running the pipeline is as follows:

*For single-end reads:*

`nextflow chipseq.nf -c ~/.nextflow/config --reads '*.fastq' --genome GRCh37 --mode single --macsconfig 'macssetup.config'`

*For paired-end reads*

`nextflow chipseq.nf -c ~/.nextflow/config --reads '*.R{1,2}.fastq' --genome GRCh37 --mode paired --macsconfig 'macssetup.config'`

##### -c

Path to the configuration file for Nextflow.

##### --reads

Location of the input FastQ files:
Note that the {1,2} parentheses are required to specify paired end data. Also, note that the file path should be in quotation marks to prevent shell glob expansion. If left unspecified, the pipeline will assume that the data is in a directory called data in the working directory.

##### --genome

The reference genome to use of the analysis, needs to be one of the genome specified in the config file. The human GRCh37 genome is set as default. Options: GRCh37, GRCm38, NCBIM37

##### --mode

single or paired

##### --macsconfig

The setup file for peak calling using MACS. 

Format:

`ChIPSampleID1,CtrlSampleID1,AnalysisID1`

`ChIPSampleID2,CtrlSampleID2,AnalysisID2`

`ChIPSampleID3,,AnalysisID3`

Note! For single-sample peaking calling without a control sample, skip the field of `CtrlSampleID`.
