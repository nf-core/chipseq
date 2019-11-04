# nf-core/chipseq: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2019-11-05

### `Added`

* [#46](https://github.com/nf-core/atacseq/issues/46) - Missing gene_bed path in igenomes config
* Update template to tools `1.7`
* Add `--trim_nextseq` parameter
* Add `CITATIONS.md` file
* Capitalised process names

### `Fixed`

* **Change all parameters from `camelCase` to `snake_case` (see [Deprecated](#Deprecated))**
* [#44](https://github.com/nf-core/atacseq/issues/44) - Output directory missing: macs2/consensus/deseq2
* [#45](https://github.com/nf-core/atacseq/issues/45) - Wrong x-axis scale for the HOMER: Peak annotation Counts tab plot?
* [#46](https://github.com/nf-core/atacseq/issues/46) - Stage blacklist file in channel properly
* [#50](https://github.com/nf-core/atacseq/issues/50) - HOMER number of peaks does not correspond to found MACS2 peaks
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
