## nf-core/chipseq v1.0dev
v1.0dev nf-core/chipseq
* Add option for building BWA index for larger ref
* Update software versions in uppmax-modules
* Add function to validate input files in MACS config file
* Remove ngsplot and move the function of plotProfile in deepTools
* Support reference genome GRCm37
* Add pre-defined genome sizes for all reference genomes to support macs2 peak calling and downstream processing
* Add blacklist files for ce10, dm3, hg38, and mm9
* Add option for special analysis setup for ATAC-seq analysis
* Documents revised accordingly.

05/02/2019 17:11:13
* Major overhaul of docs and assets in-line with nf-core/tools v1.4
* Added ability to use nf-core/configs along with associated docs
* Updated manifest scope to deal with pipeline version
* Removed NGI and SciLifeLab logos, and changed name of pipeline logo
* Added awsbatch configuration
* Put file() calls in fromFilePath()
* Removed --project param specific to UPPMAX
* Moved appropriate default params variables to nextflow.config
* Changed Picard memory specification
* Changed version number back to 1.0dev from 1.0
* Updated conda packages

Repository moved from https://github.com/SciLifeLab/NGI-ChIPseq
