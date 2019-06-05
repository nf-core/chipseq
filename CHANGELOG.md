# nf-core/chipseq: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## v1.0.0 - [2019-06-05]
Initial release of nf-core/chipseq.

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
