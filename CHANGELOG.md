# nf-core/chipseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.2.0dev - [date]

### Enhancements & fixes

- [[#428](https://github.com/nf-core/chipseq/issues/428)] - Bump MultiQC version to `1.25.1`.
- [[PR #434](https://github.com/nf-core/chipseq/pull/434)] - Prevent pipeline fails from erroneous param validation when igenomes is used.
- [[#432](https://github.com/nf-core/chipseq/issues/432)] - Fix `GFFREAD` call to have the two expected input channels.
- [[PR #444](https://github.com/nf-core/chipseq/pull/444)] - Add empty map to ch_gff so that when provided by the user `GFFREAD` works.
- Updated pipeline template to [nf-core/tools 3.2.0](https://github.com/nf-core/tools/releases/tag/3.2.0).
- [[#451](https://github.com/nf-core/chipseq/issues/451)] - Pass `map.single_read` to `SUBREAD_FEATURECOUNTS` as to correctly set parameter `-p`.

### Parameters

| Old parameter | New parameter |
| ------------- | ------------- |
|               |               |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
|            |             |             |
|            |             |             |

## [[2.1.0](https://github.com/nf-core/chipseq/releases/tag/2.1.0)] - 2024-10-07

### Credits

Special thanks to the following people for their contributions to this release:

- [Adam Talbot](https://github.com/adamrtalbot)
- [Björn Langer](https://github.com/bjlang)
- [Konrad Rokicki](https://github.com/krokicki)
- [Matthias Hörtenhuber](https://github.com/mashehu)
- [Maxime Garcia](https://github.com/maxulysse)
- [Samuel Ruiz Pérez](https://github.com/samuelruizperez)
- [Sarah Guinchard](https://github.com/g-sarah)
- [Sateesh Peri](https://github.com/sateeshperi)
- [Steffen Möller](https://github.com/smoe)

Thank you to everyone else that has contributed by reporting bugs, enhancements or in any other way, shape or form.

### Enhancements & fixes

- Updated pipeline template to [nf-core/tools 2.7.2](https://github.com/nf-core/tools/releases/tag/2.7.2)
- [[#317](https://github.com/nf-core/chipseq/issues/317)] - Added metro map
- [[#288](https://github.com/nf-core/chipseq/issues/291)] - Bump `chromap` version 2 and enable all the steps below chromap again when paired-end data is processed.
- [[#311](https://github.com/nf-core/chipseq/issues/311)] - Add back `--skip_spp` parameter which was unintentionally removed from the code.
- Install available nf-core subworkflows and refactor code accordingly
- [[#318](https://github.com/nf-core/chipseq/issues/318)] - Update `bowtie2/align` module to fix issue when downloading its singularity image.
- [[#320](https://github.com/nf-core/chipseq/issues/320)] - Fix samplesheet control column in documentation examples.
- [[#328](https://github.com/nf-core/chipseq/issues/328)] - Modify documentation to clarify that is necessary to provide the `--read_length` when `--genome` is set and `--macs_gsize` has not provided.
- Remove `enable_conda` param from local modules.
- Fix the path where `chromap` index is stored when `--save_reference` is set.
- Fix untar of `chromap` index when using `--chromap_index` param.
- [nf-core/tools#2286](https://github.com/nf-core/tools/issues/2286) - Set default container registry outside profile scope.
- [[#343](https://github.com/nf-core/chipseq/issues/343)] - Provide replicate information explicitly in samplesheet.
- Updated pipeline template to [nf-core/tools 2.10](https://github.com/nf-core/tools/releases/tag/2.10).
- [[#367](https://github.com/nf-core/chipseq/issues/367)] - Get rid of `CheckIfExists` for params paths.
- [[#370](https://github.com/nf-core/chipseq/issues/370)] - Fix stack overflow exceptions in phantompeakqualtools ([see here](https://github.com/kundajelab/phantompeakqualtools/issues/3)).
- [[#387](https://github.com/nf-core/chipseq/issues/387)] - Get rid of the `lib` folder and rearrange the pipeline accordingly.
- [[#385](https://github.com/nf-core/chipseq/issues/385)] - Fix `--save_unaligned` description in schema.
- [[PR #392](https://github.com/nf-core/chipseq/pull/392)] - Adding line numbers to warnings/errors messages in `bin/check_samplesheet.py`.
- [[#396](https://github.com/nf-core/chipseq/issues/396)] - Check that samplesheet samples IDs do only have alphanumeric characters, dots, dashes or underscores.
- [[#378](https://github.com/nf-core/chipseq/issues/378)] - Switch from macs2 to macs3.
- [[#347](https://github.com/nf-core/chipseq/issues/347)] - Add read group tag to bam files processed by bowtie2.
- [[PR #406](https://github.com/nf-core/chipseq/pull/406)] - Update metro map to show macs3 instead of macs2.
- [[#409](https://github.com/nf-core/chipseq/issues/409)] - Bulk modules and subworkflows update.
- [[PR #415](https://github.com/nf-core/chipseq/pull/415)] - Get rid of `oras` in modules.
- [[PR #420](https://github.com/nf-core/chipseq/pull/420)] - Fix bug on how replicates are set for MACS3 downstream analysis.

### Parameters

| Old parameter          | New parameter                        |
| ---------------------- | ------------------------------------ |
| `--show_hidden_params` | `--validationShowHiddenParams`       |
|                        | `--version`                          |
|                        | `--hook_url`                         |
|                        | `--multiqc_logo`                     |
|                        | `--multiqc_methods_description`      |
|                        | `--pipelines_testdata_base_path`     |
|                        | `--validationFailUnrecognisedParams` |
|                        | `--validationLenientMode`            |
| `--enable_conda`       |                                      |

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency              | Old version | New version |
| ----------------------- | ----------- | ----------- |
| `bowtie2`               | 2.4.4       | 2.5.2       |
| `bwa`                   | 0.7.17      | 0.7.18      |
| `chromap`               | 0.2.1       | 0.2.6       |
| `deeptools`             | 3.5.1       | 3.5.5       |
| `fastqc`                | 0.11.9      | 0.12.1      |
| `gffread`               | 0.12.1      | 0.12.7      |
| `macs2`                 | 2.2.7.1     |             |
| `macs3`                 |             | 3.0.1       |
| `multiqc`               | 1.13        | 1.23        |
| `picard`                | 2.27.4      | 3.2.0       |
| `samtools`              | 1.15.1      | 1.20        |
| `ucsc-bedgraphtobigwig` | 377         | 445         |
| `umi_tools`             |             | 1.1.5       |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [[2.0.0](https://github.com/nf-core/chipseq/releases/tag/2.0.0)] - 2022-10-03

### Enhancements & fixes

- Pipeline has been re-implemented in [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)
- All software containers are now exclusively obtained from [Biocontainers](https://biocontainers.pro/#/registry)
- Updated pipeline template to [nf-core/tools 2.5.1](https://github.com/nf-core/tools/releases/tag/2.5.1)
- [[#128](https://github.com/nf-core/chipseq/issues/128)] - Filter files with no peaks to avoid errors in downstream processes
- [[#220](https://github.com/nf-core/chipseq/issues/220)] - Fix `phantompeakqualtools` protection stack overflow error
- [[#233](https://github.com/nf-core/chipseq/issues/233)] - Add `chromap` to the available aligners
- Bump minimum Nextflow version from `21.04.0` -> `21.10.3`
- Added `python3` shebang to appropriate scripts in `bin/` directory
- [[#160](https://github.com/nf-core/chipseq/issues/160)] - Add `bowtie2` and `star` as available aligners, via the `--aligner` parameter
- Add `--save_unaligned` parameter (only available for `bowtie2` and `star`)
- Update `igenomes.config` to fetch whole `BWAIndex/version0.6.0/` folder
- [[228](https://github.com/nf-core/chipseq/issues/228)] - Update blacklist bed files.
- [nf-core/tools#1415](https://github.com/nf-core/tools/issues/1415) - Make `--outdir` a mandatory parameter
- [[282](https://github.com/nf-core/chipseq/issues/282)] - Fix `genome.fa` publication for IGV.
- [[280](https://github.com/nf-core/chipseq/issues/280)] - Update `macs_gsize` in `igenomes.config`, create a new `--read_length` parameter and implement the logic to calculate `--macs_gsize` when the parameter is missing
- Eliminate `if` conditions from `deseq2_qc` and `macs2_consensus` (local module and use `ext.when` instead)
- Remove `deseq2` differential binding analysis of consensus peaks.
- [[280](https://github.com/nf-core/chipseq/issues/291) - Filter paired-end files produced by `chromap` since the resulting `BAM` files can not be processed downstream.
- Add bytesize link to readme.

### Parameters

| Old parameter          | New parameter           |
| ---------------------- | ----------------------- |
| `--conda`              | `--enable_conda`        |
| `--skip_diff_analysis` | `--skip_deseq2_qc`      |
|                        | `--skip_qc`             |
|                        | `--aligner`             |
|                        | `--save_unaligned`      |
|                        | `--read_length`         |
|                        | `--multiqc_title`       |
|                        | `--gff`                 |
|                        | `--bowtie2_index`       |
|                        | `--chromap_index`       |
|                        | `--star_index`          |
|                        | `--validate_params`     |
|                        | `--show_hidden_params`  |
|                        | `--config_profile_name` |
| `--clusterOptions`     |                         |
| `--single_end`         |                         |
| `--name`               |                         |
| `--hostnames`          |                         |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency              | Old version | New version |
| ----------------------- | ----------- | ----------- |
| `samtools`              | 1.10        | 1.15.1      |
| `picard`                | 2.23.1      | 2.27.4      |
| `bamtools`              | 2.5.1       | 2.5.2       |
| `pysam`                 | 0.15.3      | 0.19.0      |
| `bedtools`              | 2.29.2      | 2.30.0      |
| `ucsc-bedgraphtobigwig` | 357         | 377         |
| `deeptools`             | 3.4.3       | 3.5.1       |
| `pigz`                  | 2.3.4       | 2.6         |
| `preseq`                | 2.0.3       | 3.1.2       |
| `multiqc`               | 1.9         | 1.13a       |
| `r-base`                | 3.6.1       | 4.0.3       |
| `r-ggplot2`             | 3.3.2       | 3.3.3       |
| `bioconductor-deseq2`   | 1.26.0      | 1.28.0      |
| `trim-galore`           | 0.6.5       | 0.6.7       |
| `r-optparse`            | -           | 1.7.1       |
| `chromap`               | -           | 0.2.1       |
| `bowtie2`               | -           | 2.4.4       |
| `star`                  | -           | 2.6.1d      |
| `r-tidyr`               | -           | -           |
| `r-lattice`             | -           | -           |
| `r-xfun`                | -           | -           |
| `bioconductor-vsn`      | -           | -           |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [[1.2.2](https://github.com/nf-core/chipseq/releases/tag/1.2.2)] - 2021-04-22

- [#206](https://github.com/nf-core/chipseq/issues/206) - Minor patch release to fix Conda environment

### `Dependencies`

- Update r-base `3.6.2` -> `3.6.3`
- Update r-xfun `0.15` -> `0.20`

## [[1.2.1](https://github.com/nf-core/chipseq/releases/tag/1.2.1)] - 2020-07-29

- [#171](https://github.com/nf-core/chipseq/issues/171) - Minor patch release to update pipeline schema

## [[1.2.0](https://github.com/nf-core/chipseq/releases/tag/1.2.0)] - 2020-07-02

### `Added`

- [#138](https://github.com/nf-core/chipseq/issues/138) - Add social preview image
- [#153](https://github.com/nf-core/chipseq/issues/153) - Add plotHeatmap
- [#159](https://github.com/nf-core/chipseq/issues/159) - expose bwa mem -T parameter
- [nf-core/atacseq#63](https://github.com/nf-core/atacseq/issues/63) - Added multicore support for Trim Galore!
- [nf-core/atacseq#75](https://github.com/nf-core/atacseq/issues/75) - Include gene annotation versions in multiqc report
- [nf-core/atacseq#76](https://github.com/nf-core/atacseq/issues/76) - featureCounts coupled to DESeq2
- [nf-core/atacseq#79](https://github.com/nf-core/atacseq/issues/79) - Parallelize DESeq2
- [nf-core/atacseq#97](https://github.com/nf-core/atacseq/issues/97) - PBC1, PBC2 from pipeline?
- [nf-core/atacseq#107](https://github.com/nf-core/atacseq/issues/107) - Add options to change MACS2 parameters
- Regenerated screenshots and added collapsible sections for output files in `docs/output.md`
- Update template to tools `1.9`
- Replace `set` with `tuple` and `file()` with `path()` in all processes
- Capitalise process names
- Parameters:
  - `--bwa_min_score` to set minimum alignment score for BWA MEM
  - `--macs_fdr` to provide FDR threshold for MACS2 peak calling
  - `--macs_pvalue` to provide p-value threshold for MACS2 peak calling
  - `--skip_peak_qc` to skip MACS2 peak QC plot generation
  - `--skip_peak_annotation` to skip annotation of MACS2 and consensus peaks with HOMER
  - `--skip_consensus_peaks` to skip consensus peak generation
  - `--deseq2_vst` to use variance stabilizing transformation (VST) instead of regularized log transformation (rlog) with DESeq2
  - `--publish_dir_mode` to customise method of publishing results to output directory [nf-core/tools#585](https://github.com/nf-core/tools/issues/585)

### `Removed`

- `--tss_bed` parameter

### `Fixed`

- [#118](https://github.com/nf-core/chipseq/issues/118) - Running on with SGE
- [#132](https://github.com/nf-core/chipseq/issues/132) - BigWig Error: sort: cannot create temporary file in '': Read-only file system
- [#154](https://github.com/nf-core/chipseq/issues/154) - computeMatrix.val.mat.gz files not zipped
- [nf-core/atacseq#71](https://github.com/nf-core/atacseq/issues/71) - consensus_peaks.mLb.clN.boolean.intersect.plot.pdf not generated
- [nf-core/atacseq#73](https://github.com/nf-core/atacseq/issues/73) - macs_annotatePeaks.mLb.clN.summary.txt file is not created
- [nf-core/atacseq#86](https://github.com/nf-core/atacseq/issues/86) - bug in the plot_homer_annotatepeaks.r script
- [nf-core/atacseq#102](https://github.com/nf-core/atacseq/issues/102) - Incorrect Group ID assigned by featurecounts_deseq2.r
- [nf-core/atacseq#109](https://github.com/nf-core/atacseq/issues/109) - Specify custom gtf but gene bed is not generated from that gtf?
- Make executables in `bin/` compatible with Python 3

### `Dependencies`

- Add bioconductor-biocparallel `1.20.0`
- Add markdown `3.2.2`
- Add pigz `2.3.4`
- Add pygments `2.6.1`
- Add pymdown-extensions `7.1`
- Add python `3.7.6`
- Add r-reshape2 `1.4.4`
- Add r-tidyr `1.1.0`
- Update bedtools `2.27.1` -> `2.29.2`
- Update bioconductor-deseq2 `1.20.0` -> `1.26.0`
- Update bioconductor-vsn `3.46.0` -> `3.54.0`
- Update deeptools `3.2.1` -> `3.4.3`
- Update fastqc `0.11.8` -> `0.11.9`
- Update gawk `4.2.1` -> `5.1.0`
- Update homer `4.9.1` -> `4.11`
- Update macs2 `2.1.2` -> `2.2.7.1`
- Update multiqc `1.7` -> `1.8`
- Update phantompeakqualtools `1.2` -> `1.2.2`
- Update picard `2.19.0` -> `2.23.1`
- Update pysam `0.15.2` -> `0.15.3`
- Update r-base `3.4.1` -> `3.6.2`
- Update r-ggplot2 `3.1.0` -> `3.3.2`
- Update r-lattice `0.20_35` -> `0.20_41`
- Update r-optparse `1.6.0` -> `1.6.6`
- Update r-pheatmap `1.0.10` -> `1.0.12`
- Update r-scales `1.0.0` -> `1.1.1`
- Update r-upsetr `1.3.3` -> `1.4.0`
- Update r-xfun `0.3` -> `0.15`
- Update samtools `1.9` -> `1.10`
- Update subread `1.6.4` -> `2.0.1`
- Update trim-galore `0.5.0` -> `0.6.5`
- Update ucsc-bedgraphtobigwig `377` -> `357`

## [[1.1.0](https://github.com/nf-core/chipseq/releases/tag/1.1.0)] - 2019-11-05

### `Added`

- [nf-core/atacseq#46](https://github.com/nf-core/atacseq/issues/46) - Missing gene_bed path in igenomes config
- Update template to tools `1.7`
- Add `--trim_nextseq` parameter
- Add `CITATIONS.md` file
- Capitalised process names

### `Fixed`

- Change all parameters from `camelCase` to `snake_case` (see [Deprecated](#Deprecated))
- [nf-core/atacseq#44](https://github.com/nf-core/atacseq/issues/44) - Output directory missing: macs2/consensus/deseq2
- [nf-core/atacseq#45](https://github.com/nf-core/atacseq/issues/45) - Wrong x-axis scale for the HOMER: Peak annotation Counts tab plot?
- [nf-core/atacseq#46](https://github.com/nf-core/atacseq/issues/46) - Stage blacklist file in channel properly
- [nf-core/atacseq#50](https://github.com/nf-core/atacseq/issues/50) - HOMER number of peaks does not correspond to found MACS2 peaks
- Fixed bug in UpSetR peak intersection plot
- Increase default resource requirements in `base.config`
- Increase process-specific requirements based on user-reported failures

### `Dependencies`

- Update Nextflow `0.32.0` -> `19.10.0`

### `Deprecated`

| Deprecated                   | Replacement               |
| ---------------------------- | ------------------------- |
| `--design`                   | `--input`                 |
| `--singleEnd`                | `--single_end`            |
| `--saveGenomeIndex`          | `--save_reference`        |
| `--skipTrimming`             | `--skip_trimming`         |
| `--saveTrimmed`              | `--save_trimmed`          |
| `--keepDups`                 | `--keep_dups`             |
| `--keepMultiMap`             | `--keep_multi_map`        |
| `--saveAlignedIntermediates` | `--save_align_intermeds`  |
| `--narrowPeak`               | `--narrow_peak`           |
| `--saveMACSPileup`           | `--save_macs_pileup`      |
| `--skipDiffAnalysis`         | `--skip_diff_analysis`    |
| `--skipFastQC`               | `--skip_fastqc`           |
| `--skipPicardMetrics`        | `--skip_picard_metrics`   |
| `--skipPreseq`               | `--skip_preseq`           |
| `--skipPlotProfile`          | `--skip_plot_profile`     |
| `--skipPlotFingerprint`      | `--skip_plot_fingerprint` |
| `--skipSpp`                  | `--skip_spp`              |
| `--skipIGV`                  | `--skip_igv`              |
| `--skipMultiQC`              | `--skip_multiqc`          |

## [[1.0.0](https://github.com/nf-core/chipseq/releases/tag/1.0.0)] - 2019-06-06

Initial release of nf-core/chipseq pipeline.

### `Added`

- Raw read QC (FastQC)
- Adapter trimming (Trim Galore!)
- Map and filter reads (BWA, picard, SAMtools, BEDTools, BAMTools, Pysam)
- Create library-size normalised bigWig tracks (BEDTools, bedGraphToBigWig)
- Alignment QC metrics (Preseq, picard)
- ChIP-seq QC metrics (deepTools, phantompeakqualtools)
- Call and annotate broad/narrow peaks (MACS2, HOMER)
- Create consensus set of peaks per antibody (BEDTools)
- Quantification and differential binding analysis (featureCounts, DESeq2)
- Collate appropriate files for genome browser visualisation (IGV)
- Collate and present various QC metrics (MultiQC, R)
