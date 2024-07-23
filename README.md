<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-chipseq_logo_dark.png">
    <img alt="nf-core/chipseq" src="docs/images/nf-core-chipseq_logo_light.png">
  </picture>
</h1>
[![GitHub Actions CI Status](https://github.com/nf-core/chipseq/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/chipseq/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/chipseq/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/chipseq/actions?query=workflow%3A%22nf-core+linting%22)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/chipseq/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.3240506-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.3240506)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/chipseq)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23chipseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/chipseq)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nfcore/chipseq** is a bioinformatics analysis pipeline used for Chromatin ImmunoPrecipitation sequencing (ChIP-seq) data.

On release, automated continuous integration tests run the pipeline on a [full-sized dataset](https://github.com/nf-core/test-datasets/tree/chipseq#full-test-dataset-origin) on the AWS cloud infrastructure. The dataset consists of FoxA1 (transcription factor) and EZH2 (histone,mark) IP experiments from _Franco et al. 2015_ ([GEO: GSE59530](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59530), [PMID: 25752574](https://pubmed.ncbi.nlm.nih.gov/25752574/)) and _Popovic et al. 2014_ ([GEO: GSE57632](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57632), [PMID: 25188243](https://pubmed.ncbi.nlm.nih.gov/25188243/)), respectively. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from running the full-sized tests can be viewed on the [nf-core website](https://nf-co.re/chipseq/results).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Online videos

A short talk about the history, current status and functionality on offer in this pipeline was given by Jose Espinosa-Carrasco ([@joseespinosa](https://github.com/joseespinosa)) on [26th July 2022](https://nf-co.re/events/2022/bytesize-chipseq) as part of the nf-core/bytesize series.

You can find numerous talks on the [nf-core events page](https://nf-co.re/events) from various topics including writing pipelines/modules in Nextflow DSL2, using nf-core tooling, running nf-core pipelines as well as more generic content like contributing to Github. Please check them out!

## Pipeline summary

![nf-core/chipseq metro map](docs/images/nf-core-chipseq_metro_map_grey.png)

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. Choice of multiple aligners
   1.([`BWA`](https://sourceforge.net/projects/bio-bwa/files/))
   2.([`Chromap`](https://github.com/haowenz/chromap))
   3.([`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
   4.([`STAR`](https://github.com/alexdobin/STAR))
4. Mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
5. Merge alignments from multiple libraries of the same sample ([`picard`](https://broadinstitute.github.io/picard/))
   1. Re-mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
   2. Filtering to remove:
      - reads mapping to blacklisted regions ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/), [`BEDTools`](https://github.com/arq5x/bedtools2/))
      - reads that are marked as duplicates ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
      - reads that are not marked as primary alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
      - reads that are unmapped ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
      - reads that map to multiple locations ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
      - reads containing > 4 mismatches ([`BAMTools`](https://github.com/pezmaster31/bamtools))
      - reads that have an insert size > 2kb ([`BAMTools`](https://github.com/pezmaster31/bamtools); _paired-end only_)
      - reads that map to different chromosomes ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); _paired-end only_)
      - reads that arent in FR orientation ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); _paired-end only_)
      - reads where only one read of the pair fails the above criteria ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); _paired-end only_)
   3. Alignment-level QC and estimation of library complexity ([`picard`](https://broadinstitute.github.io/picard/), [`Preseq`](http://smithlabresearch.org/software/preseq/))
   4. Create normalised bigWig files scaled to 1 million mapped reads ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
   5. Generate gene-body meta-profile from bigWig files ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html))
   6. Calculate genome-wide IP enrichment relative to control ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html))
   7. Calculate strand cross-correlation peak and ChIP-seq quality measures including NSC and RSC ([`phantompeakqualtools`](https://github.com/kundajelab/phantompeakqualtools))
   8. Call broad/narrow peaks ([`MACS3`](https://github.com/macs3-project/MACS))
   9. Annotate peaks relative to gene features ([`HOMER`](http://homer.ucsd.edu/homer/download.html))
   10. Create consensus peakset across all samples and create tabular file to aid in the filtering of the data ([`BEDTools`](https://github.com/arq5x/bedtools2/))
   11. Count reads in consensus peaks ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
   12. PCA and clustering ([`R`](https://www.r-project.org/), [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
6. Create IGV session file containing bigWig tracks, peaks and differential sites for data visualisation ([`IGV`](https://software.broadinstitute.org/software/igv/)).
7. Present QC for raw read, alignment, peak-calling and differential binding results ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

To run on your data, prepare a tab-separated samplesheet with your input data. Please follow the [documentation on samplesheets](https://nf-co.re/chipseq/usage#samplesheet-input) for more details. An example samplesheet for running the pipeline looks as follows:

```csv
sample,fastq_1,fastq_2,antibody,control
WT_BCATENIN_IP_REP1,BLA203A1_S27_L006_R1_001.fastq.gz,,BCATENIN,WT_INPUT_REP1
WT_BCATENIN_IP_REP2,BLA203A25_S16_L001_R1_001.fastq.gz,,BCATENIN,WT_INPUT_REP2
WT_BCATENIN_IP_REP2,BLA203A25_S16_L002_R1_001.fastq.gz,,BCATENIN,WT_INPUT_REP2
WT_BCATENIN_IP_REP2,BLA203A25_S16_L003_R1_001.fastq.gz,,BCATENIN,WT_INPUT_REP2
WT_BCATENIN_IP_REP3,BLA203A49_S40_L001_R1_001.fastq.gz,,BCATENIN,WT_INPUT_REP3
WT_INPUT_REP1,BLA203A6_S32_L006_R1_001.fastq.gz,,,
WT_INPUT_REP2,BLA203A30_S21_L001_R1_001.fastq.gz,,,
WT_INPUT_REP2,BLA203A30_S21_L002_R1_001.fastq.gz,,,
WT_INPUT_REP3,BLA203A31_S21_L003_R1_001.fastq.gz,,,
```

Now, you can run the pipeline using:

```bash
nextflow run nf-core/chipseq --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
```

See [usage docs](https://nf-co.re/chipseq/usage) for all of the available options when running the pipeline.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/chipseq/usage) and the [parameter documentation](https://nf-co.re/chipseq/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/chipseq/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/chipseq/output).

## Credits

These scripts were originally written by Chuan Wang ([@chuan-wang](https://github.com/chuan-wang)) and Phil Ewels ([@ewels](https://github.com/ewels)) for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/) at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden. The pipeline was re-implemented by Harshil Patel ([@drpatelh](https://github.com/drpatelh)) from [Seqera Labs, Spain](https://seqera.io/) and converted to Nextflow DSL2 by Jose Espinosa-Carrasco ([@JoseEspinosa](https://github.com/JoseEspinosa)) from [The Comparative Bioinformatics Group](https://www.crg.eu/en/cedric_notredame) at [The Centre for Genomic Regulation, Spain](https://www.crg.eu/).

The pipeline workflow diagram was designe by Sarah Guinchard ([@G-Sarah](https://github.com/G-Sarah)).

Many thanks to others who have helped out and contributed along the way too, including (but not limited to): [@apeltzer](https://github.com/apeltzer), [@bc2zb](https://github.com/bc2zb), [@crickbabs](https://github.com/crickbabs), [@drejom](https://github.com/drejom), [@houghtos](https://github.com/houghtos), [@KevinMenden](https://github.com/KevinMenden), [@mashehu](https://github.com/mashehu), [@pditommaso](https://github.com/pditommaso), [@Rotholandus](https://github.com/Rotholandus), [@sofiahaglund](https://github.com/sofiahaglund), [@tiagochst](https://github.com/tiagochst) and [@winni2k](https://github.com/winni2k).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#chipseq` channel](https://nfcore.slack.com/channels/chipseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/chipseq for your analysis, please cite it using the following doi: [10.5281/zenodo.3240506](https://doi.org/10.5281/zenodo.3240506)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
