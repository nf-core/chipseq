# nf-core/chipseq: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [1.2.1] - 2020-07-29

* [#171](https://github.com/nf-core/chipseq/issues/171) - Minor patch release to update pipeline schema

## [1.2.0] - 2020-07-02

### `Added`

* [#138](https://github.com/nf-core/chipseq/issues/138) - Add social preview image
* [#153](https://github.com/nf-core/chipseq/issues/153) - Add plotHeatmap
* [#159](https://github.com/nf-core/chipseq/issues/159) - expose bwa mem -T parameter
* [nf-core/atacseq#63](https://github.com/nf-core/atacseq/issues/63) - Added multicore support for Trim Galore!
* [nf-core/atacseq#75](https://github.com/nf-core/atacseq/issues/75) - Include gene annotation versions in multiqc report
* [nf-core/atacseq#76](https://github.com/nf-core/atacseq/issues/76) - featureCounts coupled to DESeq2
* [nf-core/atacseq#79](https://github.com/nf-core/atacseq/issues/79) - Parallelize DESeq2
* [nf-core/atacseq#97](https://github.com/nf-core/atacseq/issues/97) - PBC1, PBC2 from pipeline?
* [nf-core/atacseq#107](https://github.com/nf-core/atacseq/issues/107) - Add options to change MACS2 parameters
* Regenerated screenshots and added collapsible sections for output files in `docs/output.md`
* Update template to tools `1.9`
* Replace `set` with `tuple` and `file()` with `path()` in all processes
* Capitalise process names
* Parameters:
    * `--bwa_min_score` to set minimum alignment score for BWA MEM
    * `--macs_fdr` to provide FDR threshold for MACS2 peak calling
    * `--macs_pvalue` to provide p-value threshold for MACS2 peak calling
    * `--skip_peak_qc` to skip MACS2 peak QC plot generation
    * `--skip_peak_annotation` to skip annotation of MACS2 and consensus peaks with HOMER
    * `--skip_consensus_peaks` to skip consensus peak generation
    * `--deseq2_vst` to use variance stabilizing transformation (VST) instead of regularized log transformation (rlog) with DESeq2
    * `--publish_dir_mode` to customise method of publishing results to output directory [nf-core/tools#585](https://github.com/nf-core/tools/issues/585)

### `Removed`

* `--tss_bed` parameter

### `Fixed`

* [#118](https://github.com/nf-core/chipseq/issues/118) - Running on with SGE
* [#132](https://github.com/nf-core/chipseq/issues/132) - BigWig Error: sort: cannot create temporary file in '': Read-only file system
* [#154](https://github.com/nf-core/chipseq/issues/154) - computeMatrix.val.mat.gz files not zipped
* [nf-core/atacseq#71](https://github.com/nf-core/atacseq/issues/71) - consensus_peaks.mLb.clN.boolean.intersect.plot.pdf not generated
* [nf-core/atacseq#73](https://github.com/nf-core/atacseq/issues/73) - macs_annotatePeaks.mLb.clN.summary.txt file is not created
* [nf-core/atacseq#86](https://github.com/nf-core/atacseq/issues/86) - bug in the plot_homer_annotatepeaks.r script
* [nf-core/atacseq#102](https://github.com/nf-core/atacseq/issues/102) - Incorrect Group ID assigned by featurecounts_deseq2.r
* [nf-core/atacseq#109](https://github.com/nf-core/atacseq/issues/109) - Specify custom gtf but gene bed is not generated from that gtf?
* Make executables in `bin/` compatible with Python 3

### `Dependencies`

* Add bioconductor-biocparallel `1.20.0`
* Add markdown `3.2.2`
* Add pigz `2.3.4`
* Add pygments `2.6.1`
* Add pymdown-extensions `7.1`
* Add python `3.7.6`
* Add r-reshape2 `1.4.4`
* Add r-tidyr `1.1.0`
* Update bedtools `2.27.1` -> `2.29.2`
* Update bioconductor-deseq2 `1.20.0` -> `1.26.0`
* Update bioconductor-vsn `3.46.0` -> `3.54.0`
* Update deeptools `3.2.1` -> `3.4.3`
* Update fastqc `0.11.8` -> `0.11.9`
* Update gawk `4.2.1` -> `5.1.0`
* Update homer `4.9.1` -> `4.11`
* Update macs2 `2.1.2` -> `2.2.7.1`
* Update multiqc `1.7` -> `1.8`
* Update phantompeakqualtools `1.2` -> `1.2.2`
* Update picard `2.19.0` -> `2.23.1`
* Update pysam `0.15.2` -> `0.15.3`
* Update r-base `3.4.1` -> `3.6.2`
* Update r-ggplot2 `3.1.0` -> `3.3.2`
* Update r-lattice `0.20_35` -> `0.20_41`
* Update r-optparse `1.6.0` -> `1.6.6`
* Update r-pheatmap `1.0.10` -> `1.0.12`
* Update r-scales `1.0.0` -> `1.1.1`
* Update r-upsetr `1.3.3` -> `1.4.0`
* Update r-xfun `0.3` -> `0.15`
* Update samtools `1.9` -> `1.10`
* Update subread `1.6.4` -> `2.0.1`
* Update trim-galore `0.5.0` -> `0.6.5`
* Update ucsc-bedgraphtobigwig `377` -> `357`

## [1.1.0] - 2019-11-05

### `Added`

* [nf-core/atacseq#46](https://github.com/nf-core/atacseq/issues/46) - Missing gene_bed path in igenomes config
* Update template to tools `1.7`
* Add `--trim_nextseq` parameter
* Add `CITATIONS.md` file
* Capitalised process names

### `Fixed`

* **Change all parameters from `camelCase` to `snake_case` (see [Deprecated](#Deprecated))**
* [nf-core/atacseq#44](https://github.com/nf-core/atacseq/issues/44) - Output directory missing: macs2/consensus/deseq2
* [nf-core/atacseq#45](https://github.com/nf-core/atacseq/issues/45) - Wrong x-axis scale for the HOMER: Peak annotation Counts tab plot?
* [nf-core/atacseq#46](https://github.com/nf-core/atacseq/issues/46) - Stage blacklist file in channel properly
* [nf-core/atacseq#50](https://github.com/nf-core/atacseq/issues/50) - HOMER number of peaks does not correspond to found MACS2 peaks
* Fixed bug in UpSetR peak intersection plot
* Increase default resource requirements in `base.config`
* Increase process-specific requirements based on user-reported failures

### `Dependencies`

* Update Nextflow `0.32.0` -> `19.10.0`

### `Deprecated`

| Deprecated                   | Replacement               |
|------------------------------|---------------------------|
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

## [1.0.0] - 2019-06-06

Initial release of nf-core/chipseq pipeline.

### `Added`

* Raw read QC (FastQC)
* Adapter trimming (Trim Galore!)
* Map and filter reads (BWA, picard, SAMtools, BEDTools, BAMTools, Pysam)
* Create library-size normalised bigWig tracks (BEDTools, bedGraphToBigWig)
* Alignment QC metrics (Preseq, picard)
* ChIP-seq QC metrics (deepTools, phantompeakqualtools)
* Call and annotate broad/narrow peaks (MACS2, HOMER)
* Create consensus set of peaks per antibody (BEDTools)
* Quantification and differential binding analysis (featureCounts, DESeq2)
* Collate appropriate files for genome browser visualisation (IGV)
* Collate and present various QC metrics (MultiQC, R)
