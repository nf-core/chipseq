#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                    C H I P - S E Q   B E S T   P R A C T I C E
========================================================================================
 ChIP-seq Analysis Pipeline. Started May 2016.
 @Authors
 Chuan Wang <chuan.wang@scilifelab.se>
 Phil Ewels <phil.ewels@scilifelab.se>
 Release 2016-09-23
----------------------------------------------------------------------------------------
 Basic command:
 $ nextflow main.nf

 Pipeline variables can be configured with the following command line options:
 --genome [ID]
 --index [path to bwa index] (alternative to --genome)
 --reads [path to input files]
 --macsconfig [path to config file for MACS: line format: chip_sample_id,ctrl_sample_id,analysis_id]

 For example:
 $ nextflow main.nf -c ~/.nextflow/config --reads '*.fastq.gz' --macsconfig 'macssetup.config'
 $ nextflow main.nf -c ~/.nextflow/config --reads '*.R{1,2}.fastq.gz' --macsconfig 'macssetup.config'
---------------------------------------------------------------------------------------
 The pipeline can determine whether the input data is single or paired end. This relies on
 specifying the input files correctly. For paired en data us the example above, i.e.
 'sample_*_{1,2}.fastq.gz'. Without the glob {1,2} (or similiar) the data will be treated
 as single end.
----------------------------------------------------------------------------------------
 Pipeline overview:
 - FastQC
 - TrimGalore
 - bwa
 - picard
 - phantompeakqualtools
 - deeptools
 - ngsplot
 - macs2
 - MultiQC
----------------------------------------------------------------------------------------
*/



/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 1.2

// Configurable variables
params.project = false
params.genome = false
params.index = params.genomes[ params.genome ].bwa
params.name = "NGI ChIP-seq Best Practice"
params.reads = "data/*{1,2}*.fastq.gz"
params.macsconfig = "data/macsconfig"
params.extendReadsLen = 100
params.outdir = './results'

// Validate inputs
index = file(params.index)
macsconfig = file(params.macsconfig)
if( !index.exists() ) exit 1, "Missing BWA index: '$index'"
if( !macsconfig.exists() ) exit 1, "Missing MACS config: '$macsconfig'. Specify path with --macsconfig"
if( workflow.profile == 'standard' && !params.project ) exit 1, "No UPPMAX project ID found! Use --project"

// R library locations
params.rlocation = false
if (params.rlocation){
    nxtflow_libs = file(params.rlocation)
    nxtflow_libs.mkdirs()
}

params.rlocation = "$HOME/R/nxtflow_libs/"
nxtflow_libs = file(params.rlocation)
nxtflow_libs.mkdirs()

// Variable initialisation
def single

// Custom trimming options
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

log.info "===================================="
log.info " ChIP-seq: v${version}"
log.info "===================================="
log.info "Reads          : ${params.reads}"
log.info "Genome         : ${params.genome}"
log.info "BWA Index      : ${params.index}"
log.info "MACS Config    : ${params.macsconfig}"
log.info "Extend Reads   : ${params.extendReadsLen} bp"
log.info "Current home   : $HOME"
log.info "Current user   : $USER"
log.info "Current path   : $PWD"
log.info "R libraries    : ${params.rlocation}"
log.info "Script dir     : $baseDir"
log.info "Working dir    : $workDir"
log.info "Output dir     : ${params.outdir}"
if( params.clip_r1 > 0) log.info "Trim R1        : ${params.clip_r1}"
if( params.clip_r2 > 0) log.info "Trim R2        : ${params.clip_r2}"
if( params.three_prime_clip_r1 > 0) log.info "Trim 3' R1     : ${params.three_prime_clip_r1}"
if( params.three_prime_clip_r2 > 0) log.info "Trim 3' R2     : ${params.three_prime_clip_r2}"
log.info "Config Profile : ${workflow.profile}"
log.info "===================================="

/*
 * Create a channel for input read files
 */
Channel
    .fromFilePairs( params.reads, size: -1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { raw_reads_fastqc; raw_reads_trimgalore }


/*
 * Create a channel for macs config file
 */
Channel
    .from(macsconfig.readLines())
    .map { line ->
        list = line.split(',')
        chip_sample_id = list[0]
        ctrl_sample_id = list[1]
        analysis_id = list[2]
        [ chip_sample_id, ctrl_sample_id, analysis_id ]
    }
    .into{ macs_para }



/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}


/*
 * STEP 2 - Trim Galore!
 */
process trim_galore {
    tag "$name"
    publishDir "${params.outdir}/trimgalore", mode: 'copy'

    input:
    set val(name), file(reads) from raw_reads_trimgalore

    output:
    file '*.fq.gz' into trimmed_reads
    file '*trimming_report.txt' into trimgalore_results
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

    script:
    single = reads instanceof Path
    c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
    if(single) {
        """
        trim_galore --fastqc --gzip $c_r1 $tpc_r1 $reads
        """
    } else {
        """
        trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
        """
    }
}


/*
 * STEP 3.1 - align with bwa
 */
process bwa {
    tag "$reads"
    publishDir "${params.outdir}/bwa", mode: 'copy'

    input:
    file reads from trimmed_reads

    output:
    file '*.bam' into bwa_bam
    stdout into bwa_logs

    script:
    prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    set -o pipefail
    bwa mem -M $index $reads | samtools view -bT $index - > ${prefix}.bam
    """
}


/*
 * STEP 3.2 - post-alignment processing
 */

process samtools {
    tag "$bwa_bam"
    publishDir "${params.outdir}/bwa", mode: 'copy'

    input:
    file bwa_bam

    output:
    file '*.sorted.bam' into bam_picard
    file '*.sorted.bam.bai' into bwa_bai
    file '*.sorted.bed' into bed_total

    script:
    prefix = bwa_bam[0].toString() - ~/(\.bam)?$/
    """
    set -o pipefail
    samtools sort ${prefix}.bam ${prefix}.sorted
    samtools index ${prefix}.sorted.bam
    bedtools bamtobed -i ${prefix}.sorted.bam | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 > ${prefix}.sorted.bed
    """
}


/*
 * STEP 4 Picard
 */

process picard {
    tag "$bam_picard"
    publishDir "${params.outdir}/picard", mode: 'copy'

    input:
    file bam_picard

    output:
    file '*.dedup.sorted.bam' into bam_dedup_spp, bam_dedup_ngsplotconfig, bam_dedup_ngsplot, bam_dedup_deepTools, bam_dedup_macs
    file '*.dedup.sorted.bam.bai' into bai_dedup_deepTools, bai_dedup_ngsplot, bai_dedup_macs
    file '*.dedup.sorted.bed' into bed_dedup
    file '*.picardDupMetrics.txt' into picard_reports

    script:
    prefix = bam_picard[0].toString() - ~/(\.sorted)?(\.bam)?$/
    """
    set -o pipefail
    java -Xmx2g -jar \$PICARD_HOME/picard.jar MarkDuplicates \\
        INPUT=$bam_picard \\
        OUTPUT=${prefix}.dedup.bam \\
        ASSUME_SORTED=true \\
        REMOVE_DUPLICATES=true \\
        METRICS_FILE=${prefix}.picardDupMetrics.txt \\
        VALIDATION_STRINGENCY=LENIENT \\
        PROGRAM_RECORD_ID='null'

    samtools sort ${prefix}.dedup.bam ${prefix}.dedup.sorted
    samtools index ${prefix}.dedup.sorted.bam
    bedtools bamtobed -i ${prefix}.dedup.sorted.bam | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 > ${prefix}.dedup.sorted.bed
    """
}


/*
 * STEP 5 Read_count_statistics
 */

process countstat {
    publishDir "${params.outdir}/countstat", mode: 'copy'

    input:
    file input from bed_total.mix(bed_dedup).toSortedList()

    output:
    file 'read_count_statistics.txt' into countstat_results

    script:
    """
    countstat.pl $input
    """
}


/*
 * STEP 6.1 Phantompeakqualtools
 */

process phantompeakqualtools {
    tag "$bam_dedup_spp"
    publishDir "${params.outdir}/phantompeakqualtools", mode: 'copy'

    input:
    file bam_dedup_spp

    output:
    file '*.pdf' into spp_results
    file '*.spp.out' into spp_out, spp_out_mqc

    script:
    prefix = bam_dedup_spp[0].toString() - ~/(\.dedup)?(\.sorted)?(\.bam)?$/
    """
    run_spp.R -c=$bam_dedup_spp -savp -out=${prefix}.spp.out
    """
}


/*
 * STEP 6.2 Combine_spp_out
 */

process combinesppout {
    publishDir "${params.outdir}/phantompeakqualtools", mode: 'copy'

    input:
    file spp_out_list from spp_out.toSortedList()

    output:
    file 'cross_correlation.txt' into cross_correlation

    script:
    """
    cat $spp_out_list > cross_correlation.txt
    """
}


/*
 * STEP 6.3 Calculate_NSC_RSC
 */

process calculateNSCRSC {
    publishDir "${params.outdir}/phantompeakqualtools", mode: 'copy'

    input:
    file cross_correlation

    output:
    file 'crosscorrelation_processed.txt' into calculateNSCRSC_results

    script:
    """
    calculateNSCRSC.r $cross_correlation
    """
}


/*
 * STEP 7 deepTools
 */

process deepTools {
    publishDir "${params.outdir}/deepTools", mode: 'copy'

    input:
    file bam from bam_dedup_deepTools.toSortedList()
    file bai from bai_dedup_deepTools.toSortedList()

    output:
    file 'multiBamSummary.npz' into deepTools_bamsummary
    file '*.{pdf,png}' into deepTools_results

    script:
    """
    plotFingerprint \\
        -b $bam \\
        --plotFile fingerprints.pdf \\
        --extendReads=${params.extendReadsLen} \\
        --skipZeros \\
        --ignoreDuplicates \\
        --numberOfSamples 50000 \\
        --binSize=500 \\
        --plotFileFormat=pdf \\
        --plotTitle="Fingerprints"

    multiBamSummary \\
        bins \\
        --binSize=10000 \\
        --bamfiles $bam \\
        -out multiBamSummary.npz \\
        --extendReads=${params.extendReadsLen} \\
        --ignoreDuplicates \\
        --centerReads

    plotCorrelation \\
        -in multiBamSummary.npz \\
        -o scatterplot_PearsonCorr_multiBamSummary.png \\
        --corMethod pearson \\
        --skipZeros \\
        --removeOutliers \\
        --plotTitle "Pearson Correlation of Read Counts" \\
        --whatToPlot scatterplot \\

    plotCorrelation \\
        -in multiBamSummary.npz \\
        -o heatmap_SpearmanCorr_multiBamSummary.png \\
        --corMethod spearman \\
        --skipZeros \\
        --plotTitle "Spearman Correlation of Read Counts" \\
        --whatToPlot heatmap \\
        --colorMap RdYlBu \\
        --plotNumbers \\
    """
}


/*
 * STEP 8.1 Generate config file for ngsplot
 */

process ngs_config_generate {
    publishDir "${params.outdir}/ngsplot", mode: 'copy'

    input:
    file input_files from bam_dedup_ngsplotconfig.flatten().toSortedList()

    output:
    file 'ngsplot_config' into ngsplot_config

    script:
    """
    ngs_config_generate.r $input_files
    """
}


/*
 * STEP 8.2 Ngsplot
 */

process ngsplot {
    publishDir "${params.outdir}/ngsplot", mode: 'copy'

    input:
    file input_bam_files from bam_dedup_ngsplot.flatten().toSortedList()
    file input_bai_files from bai_dedup_ngsplot.flatten().toSortedList()
    file ngsplot_config from ngsplot_config

    output:
    file '*.pdf' into ngsplot_results

    script:
    def REF
    if (params.genome == 'GRCh37'){ REF = 'hg19' }
    else if (params.genome == 'GRCm38'){ REF = 'mm10' }
    else { error "No reference / reference not supported available for ngsplot! >${params.genome}<" }
    """
    ngs.plot.r \\
        -G $REF \\
        -R genebody \\
        -C $ngsplot_config \\
        -O Genebody \\
        -D ensembl \\
        -FL 300

    ngs.plot.r \\
        -G $REF \\
        -R tss \\
        -C $ngsplot_config \\
        -O TSS \\
        -FL 300
    """
}


/*
 * STEP 9 MACS
 */

process macs {
    publishDir "${params.outdir}/macs", mode: 'copy'

    input:
    file bam_for_macs from bam_dedup_macs.flatten().toSortedList()
    file bai_for_macs from bai_dedup_macs.flatten().toSortedList()
    set chip_sample_id, ctrl_sample_id, analysis_id from macs_para

    output:
    file '*.{bed,xls,r,narrowPeak}' into macs_results

    script:
    def REF
    if (params.genome == 'GRCh37'){ REF = 'hs' }
    else if (params.genome == 'GRCm38'){ REF = 'mm' }
    else { error "No reference / reference not supported available for MACS! >${params.genome}<" }
    def ctrl = ctrl_sample_id == '' ? '' : "-c ${ctrl_sample_id}.dedup.sorted.bam"
    """
    macs2 callpeak \\
        -t ${chip_sample_id}.dedup.sorted.bam \\
        $ctrl \\
        -f BAM \\
        -g $REF \\
        -n $analysis_id \\
        -q 0.01
    """
}


/*
 * STEP 10 MultiQC
 */

process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file ('fastqc/*') from fastqc_results.flatten().toList()
    file ('trimgalore/*') from trimgalore_results.flatten().toList()
    file ('bwa/*') from bwa_logs.flatten().toList()
    file ('picard/*') from picard_reports.flatten().toList()
    file ('phantompeakqualtools/*') from spp_out_mqc.flatten().toList()

    output:
    file '*multiqc_report.html'
    file '*multiqc_data'

    script:
    """
    multiqc -f .
    """
}
