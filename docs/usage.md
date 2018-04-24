# nf-core/ChIPseq Usage

## Running the pipeline
The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/ChIPseq --reads '*_R{1,2}.fastq.gz' --macsconfig 'macssetup.config'
```

Note that the pipeline will create files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/ChIPseq
```

### Reproducibility
It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/ChIPseq releases page](https://github.com/nf-core/ChIPseq/releases) and find the latest version number - numeric only (eg. `1.4`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.4`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when youook back in the future.

## Input Data

### `--reads`
Location of the input FastQ files:

```bash
 --reads 'path/to/data/sample_*_{1,2}.fastq'
```

**NB: Must be enclosed in quotes!**

Note that the `{1,2}` parentheses are required to specify paired end data. The file path should be in quotation marks to prevent shell glob expansion.

If left unspecified, the pipeline will assume that the data is in a directory called `data` in the working directory (`data/*{1,2}.fastq.gz`).

### `--singleEnd`
By default, the pipeline expects paired-end data. If you have single-end data, specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example: `--singleEnd --reads '*.fastq'`

It is not possible to run a mixture of single-end and paired-end files in one run.

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
    * Typically the file **base name**, shared between both reads _(no file type extension)_
    * _eg._ for `chip_sample_1_R1.fastq.gz` and `chip_sample_1_R2.fastq.gz`, enter `chip_sample_1`
2. Column 2: Control sample name
    * Typically the input file **base name**, shared between both reads _(no file type extension)_
    * _eg._ for `chip_input_R1.fastq.gz` and `chip_input_R2.fastq.gz`, enter `chip_input`
3. Column 3: Analysis ID
    * The analysis ID. Used for the output directory name.
    * Should be unique for each sample / line in the file.

For single-sample peaking calling without a control sample, leave the second column blank
(control sample name).

### `--broad`
Run MACS with the `--broad` flag. With this flag on, MACS will try to composite broad regions in BED12 ( a gene-model-like format ) by putting nearby highly enriched regions into a broad region with loose cutoff. The broad region is controlled by the default qvalue cutoff 0.1.

### `--saturation`
Run saturation analysis by sub-sampling the sequence reads of ChIP sample from 10% to 100% with 10% interval, then calling ChIP-seq peaks with MACS. For one test there will be 10 sets of ChIP-seq peaks. A summary file (.CSV) will be provided with the numbers of ChIP-seq peaks.

### `--blacklist_filtering`
Specifying this flag instructs the pipeline to use bundled ENCODE blacklist regions to filter out known blacklisted regions in the called ChIP-seq peaks. Please note that this is only supported when --genome is set to `GRCh37` or `GRCm38`.

### `--blacklist`
If you prefer, you can specify the full path to the blacklist regions (should be in .BED format) which will be filtered out from the called ChIP-seq peaks. Please note that `--blacklist_filtering` is required for using this option.
```bash
--blacklist_filtering --blacklist '[path to blacklisted regions]'
```

### `--extendReadsLen`
Number of base pairs to extend the reads for the deepTools analysis. This number should be
based on the length of your reads and the expected fragment length.

Default: `100`

## Reference Genomes

### `--genome`
The reference genome to use for the analysis, needs to be one of the genome specified in the config file. This is `False` by default and needs to be specified (unless index files are supplied, see below).

_Currently only the human and mouse genomes are fully supported by this pipeline. For other genomes `MACS` and `NGSplot` will not run!_

* Human
  * `--genome GRCh37`
* Mouse
  * `--genome GRCm38`

See [`conf/uppmax.config`](conf/uppmax.config) for a list of the reference genomes and their keys.

If you're not running on UPPMAX (the default profile), you can create your own config file with paths to your reference genomes. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to add this.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh37' {
        bwa = '<path to the bwa index folder>'
        fasta = '<path to the fasta file>' // used if bwa index not given
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

### `--bwa_index`
If you prefer, you can specify the full path to your reference genome when you run the pipeline:
```bash
--bwa_index '[path to BWA index]'
```

### `--fasta`
If you don't have a BWA index available, you can pass a FASTA file to the pipeline and a BWA index
will be generated for you. Combine with `--saveReference` to save for future runs.
```bash
--fasta '[path to FASTA file]'
```

### `--gtf`
The full path to GTF file can be specified for annotating peaks. Note that the GTF file should be in the Ensembl format.
```bash
--gtf '[path to GTF file]'
```

### `--allow_multi_align`
Specifying `--allow_multi_align` will turn off the filtering of secondary alignments and unmapped reads from the BWA output file. Without this option, only primary alignments will be retained.

### `--saveReference`
Supply this parameter to save any generated reference genome files to your results folder. These can then be used for future pipeline runs, reducing processing times.

## Adapter Trimming
The pipeline accepts a number of parameters to change how the trimming is done, according to your data type.
You can specify custom trimming parameters as follows:

* `--clip_r1 <NUMBER>`
  * Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads).
* `--clip_r2 <NUMBER>`
  * Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only).
* `--three_prime_clip_r1 <NUMBER>`
  * Instructs Trim Galore to remove bp from the 3' end of read 1 _AFTER_ adapter/quality trimming has been
* `--three_prime_clip_r2 <NUMBER>`
  * Instructs Trim Galore to re move bp from the 3' end of read 2 _AFTER_ adapter/quality trimming has been performed.

### `--notrim`
Specifying `--notrim` will skip the adapter trimming step. Use this if your input FastQ files have already been trimmed outside of the workflow or if you're very confident that there is no adapter contamination in your data.

### `--saveTrimmed`
By default, trimmed FastQ files will not be saved to the results directory. Specify this flag (or set to true in your config file) to copy these files to the results directory when complete.

### `--saveAlignedIntermediates`
By default, intermediate BAM files will not be saved. The final BAM files created
after the Picard MarkDuplicates step are always saved. Set to true to also copy out BAM
files from BWA and sorting steps.

## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits on UPPMAX with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples.

## Other command line parameters
### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command). Useful if using different UPPMAX
projects or different sets of reference genomes.

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override defaults. For example, we run on UPPMAX but don't want to use the MultiQC
environment module as is the default. So we specify a config file using `-c` that contains the following:

```nextflow
process.$multiqc.module = []
```

### `--clusterOptions`
Submit arbitrary SLURM options (UPPMAX profile only). For instance, you could use `--clusterOptions '-p devcore'`
to run on the development node (though won't work with default process time requests).

### `--multiqc_config`
If you would like to supply a custom config file to MultiQC, you can specify a path with `--multiqc_config`. This is used instead of the config file specific to the pipeline.
