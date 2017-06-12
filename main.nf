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
----------------------------------------------------------------------------------------
 Basic command:
 $ nextflow main.nf

 Pipeline variables can be configured with the following command line options:
 --genome [ID]
 OR --index [path to bwa index] (alternative to --genome)
 OR --fasta [path to fasta reference] (builds bwa index)
 --reads [path to input files]
 --macsconfig [path to config file for MACS: line format: chip_sample_id,ctrl_sample_id,analysis_id]

 For example:
 $ nextflow run SciLifeLab/NGI-ChIPseq --reads '*.fastq.gz' --macsconfig 'macssetup.config'
 $ nextflow run SciLifeLab/NGI-ChIPseq --reads '*.R{1,2}.fastq.gz' --macsconfig 'macssetup.config'
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
params.genomes = []
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
params.name = "NGI ChIP-seq Best Practice"
params.reads = "data/*{1,2}*.fastq.gz"
params.macsconfig = "data/macsconfig"
params.extendReadsLen = 100
params.saveReference = false
params.saveTrimmed = false
params.saveAlignedIntermediates = false
params.outdir = './results'

// Validate inputs
macsconfig = file(params.macsconfig)
if( !macsconfig.exists() ) exit 1, "Missing MACS config: '$macsconfig'. Specify path with --macsconfig"
if( params.bwa_index ){
    bwa_index = Channel
        .fromPath(params.bwa_index)
        .ifEmpty { exit 1, "BWA index not found: ${params.bwa_index}" }
} else if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
} else {
    exit 1, "No reference genome specified!"
}
if( workflow.profile == 'standard' && !params.project ) exit 1, "No UPPMAX project ID found! Use --project"

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
if(params.bwa_index)  log.info "BWA Index      : ${params.bwa_index}"
else if(params.fasta) log.info "Fasta Ref      : ${params.fasta}"
log.info "MACS Config    : ${params.macsconfig}"
log.info "Extend Reads   : ${params.extendReadsLen} bp"
log.info "Current home   : $HOME"
log.info "Current user   : $USER"
log.info "Current path   : $PWD"
log.info "Script dir     : $baseDir"
log.info "Working dir    : $workDir"
log.info "Output dir     : ${params.outdir}"
log.info "Save Reference : ${params.saveReference}"
log.info "Save Trimmed   : ${params.saveTrimmed}"
log.info "Save Intermeds : ${params.saveAlignedIntermediates}"
if( params.clip_r1 > 0) log.info "Trim R1        : ${params.clip_r1}"
if( params.clip_r2 > 0) log.info "Trim R2        : ${params.clip_r2}"
if( params.three_prime_clip_r1 > 0) log.info "Trim 3' R1     : ${params.three_prime_clip_r1}"
if( params.three_prime_clip_r2 > 0) log.info "Trim 3' R2     : ${params.three_prime_clip_r2}"
log.info "Config Profile : ${workflow.profile}"
if(workflow.commitId){ log.info "Pipeline Commit: ${workflow.commitId}" }
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
 * Reference to use for MACS and ngs.plot.r
 */
def REF_macs = false
def REF_ngsplot = false
if (params.genome == 'GRCh37'){ REF_macs = 'hs'; REF_ngsplot = 'hg19' }
else if (params.genome == 'GRCm38'){ REF_macs = 'mm'; REF_ngsplot = 'mm10' }
else if (params.genome == false){
    log.warn "No reference supplied for MACS / ngs_plot. Use '--genome GRCh37' or '--genome GRCm38' to run MACS and ngs_plot."
} else {
    log.warn "Reference '${params.genome}' not supported by MACS / ngs_plot (only GRCh37 and GRCm38)."
}


/*
 * PREPROCESSING - Build BWA index
 */
if(!params.bwa_index && fasta){
    process makeBWAindex {
        tag fasta
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta

        output:
        file "BWAIndex" into bwa_index

        script:
        """
        mkdir BWAIndex
        bwa index -a bwtsw $fasta
        """
    }
}



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
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else params.saveTrimmed ? filename : null
        }

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
    tag "$prefix"
    publishDir "${params.outdir}/bwa", mode: 'copy', saveAs: {filename -> params.saveAlignedIntermediates ? filename : null }

    input:
    file reads from trimmed_reads
    file index from bwa_index.first()

    output:
    file '*.bam' into bwa_bam
    stdout into bwa_logs

    script:
    prefix = reads[0].toString() - ~/(_1)?(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    set -o pipefail
    bwa mem -M ${index}/genome.fa $reads | samtools view -bT $index - > ${prefix}.bam
    """
}


/*
 * STEP 3.2 - post-alignment processing
 */

process samtools {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/bwa", mode: 'copy', saveAs: {filename -> params.saveAlignedIntermediates ? filename : null }

    input:
    file bam from bwa_bam

    output:
    file '*.sorted.bam' into bam_picard
    file '*.sorted.bam.bai' into bwa_bai
    file '*.sorted.bed' into bed_total

    script:
    """
    set -o pipefail
    samtools sort $bam -o ${bam.baseName}.sorted.bam
    samtools index ${bam.baseName}.sorted.bam
    bedtools bamtobed -i ${bam.baseName}.sorted.bam | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 > ${bam.baseName}.sorted.bed
    """
}


/*
 * STEP 4 Picard
 */

process picard {
    tag "$prefix"
    publishDir "${params.outdir}/picard", mode: 'copy'

    input:
    file bam from bam_picard

    output:
    file '*.dedup.sorted.bam' into bam_dedup_spp, bam_dedup_ngsplot, bam_dedup_deepTools, bam_dedup_macs
    file '*.dedup.sorted.bam.bai' into bai_dedup_deepTools, bai_dedup_ngsplot, bai_dedup_macs
    file '*.dedup.sorted.bed' into bed_dedup
    file '*.picardDupMetrics.txt' into picard_reports

    script:
    prefix = bam[0].toString() - ~/(\.sorted)?(\.bam)?$/
    if( task.memory == null ){
        log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
        avail_mem = 3000
    } else {
        avail_mem = task.memory.toMega()
        if( avail_mem <= 0){
            avail_mem = 3000
            log.info "[Picard MarkDuplicates] Available memory 0 - defaulting to 3GB. Specify process memory requirements to change this."
        } else if( avail_mem < 250){
            avail_mem = 250
            log.info "[Picard MarkDuplicates] Available memory under 250MB - defaulting to 250MB. Specify process memory requirements to change this."
        }
    }
    """
    set -o pipefail
    java -Xmx${avail_mem}m -jar \$PICARD_HOME/picard.jar MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${prefix}.dedup.bam \\
        ASSUME_SORTED=true \\
        REMOVE_DUPLICATES=true \\
        METRICS_FILE=${prefix}.picardDupMetrics.txt \\
        VALIDATION_STRINGENCY=LENIENT \\
        PROGRAM_RECORD_ID='null'

    samtools sort ${prefix}.dedup.bam -o ${prefix}.dedup.sorted.bam
    samtools index ${prefix}.dedup.sorted.bam
    bedtools bamtobed -i ${prefix}.dedup.sorted.bam | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 > ${prefix}.dedup.sorted.bed
    """
}


/*
 * STEP 5 Read_count_statistics
 */

process countstat {
    tag "${input[0].baseName}"
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
    tag "$prefix"
    publishDir "${params.outdir}/phantompeakqualtools", pattern: '*.pdf', mode: 'copy'

    input:
    file bam from bam_dedup_spp

    output:
    file '*.pdf' into spp_results
    file '*.spp.out' into spp_out, spp_out_mqc

    script:
    prefix = bam[0].toString() - ~/(\.dedup)?(\.sorted)?(\.bam)?$/
    """
    script_path=\$(which run_spp.R)
    Rscript \${script_path} -c="$bam" -savp -out="${prefix}.spp.out"
    """
}


/*
 * STEP 6.2 Combine and calculate NSC & RSC
 */

process calculateNSCRSC {
    tag "${spp_out_list[0].baseName}"
    publishDir "${params.outdir}/phantompeakqualtools", mode: 'copy'

    input:
    file spp_out_list from spp_out.collect()

    output:
    file 'cross_correlation_processed.txt' into calculateNSCRSC_results

    script:
    """
    cat $spp_out_list > cross_correlation.txt
    calculateNSCRSC.r cross_correlation.txt
    """
}


/*
 * STEP 7 deepTools
 */

process deepTools {
    tag "${bam[0].baseName}"
    publishDir "${params.outdir}/deepTools", mode: 'copy'

    input:
    file bam from bam_dedup_deepTools.collect()
    file bai from bai_dedup_deepTools.collect()

    output:
    file '*.{pdf,png,npz}' into deepTools_results

    script:
    if(bam instanceof Path){
        log.warn("Only 1 BAM file - skipping multiBam deepTool steps")
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
        """
    } else {
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
            --whatToPlot scatterplot

        plotCorrelation \\
            -in multiBamSummary.npz \\
            -o heatmap_SpearmanCorr_multiBamSummary.png \\
            --corMethod spearman \\
            --skipZeros \\
            --plotTitle "Spearman Correlation of Read Counts" \\
            --whatToPlot heatmap \\
            --colorMap RdYlBu \\
            --plotNumbers
        """
    }
}


/*
 * STEP 8 Ngsplot
 */

process ngsplot {
    tag "${input_bam_files[0].baseName}"
    publishDir "${params.outdir}/ngsplot", mode: 'copy'

    input:
    file input_bam_files from bam_dedup_ngsplot.collect()
    file input_bai_files from bai_dedup_ngsplot.collect()

    output:
    file '*.pdf' into ngsplot_results

    when: REF_ngsplot

    script:
    """
    ngs_config_generate.r $input_bam_files

    ngs.plot.r \\
        -G $REF_ngsplot \\
        -R genebody \\
        -C ngsplot_config \\
        -O Genebody \\
        -D ensembl \\
        -FL 300

    ngs.plot.r \\
        -G $REF_ngsplot \\
        -R tss \\
        -C ngsplot_config \\
        -O TSS \\
        -FL 300
    """
}


/*
 * STEP 9 MACS
 */

process macs {
    tag "${bam_for_macs[0].baseName}"
    publishDir "${params.outdir}/macs", mode: 'copy'

    input:
    file bam_for_macs from bam_dedup_macs.collect()
    file bai_for_macs from bai_dedup_macs.collect()
    set chip_sample_id, ctrl_sample_id, analysis_id from macs_para

    output:
    file '*.{bed,xls,r,narrowPeak}' into macs_results

    when: REF_macs

    script:
    def ctrl = ctrl_sample_id == '' ? '' : "-c ${ctrl_sample_id}.dedup.sorted.bam"
    """
    macs2 callpeak \\
        -t ${chip_sample_id}.dedup.sorted.bam \\
        $ctrl \\
        -f BAM \\
        -g $REF_macs \\
        -n $analysis_id \\
        -q 0.01
    """
}


/*
 * STEP 10 MultiQC
 */

process multiqc {
    tag "$prefix"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file (fastqc:'fastqc/*') from fastqc_results.collect()
    file ('trimgalore/*') from trimgalore_results.collect()
    file ('bwa/*') from bwa_logs.collect()
    file ('picard/*') from picard_reports.collect()
    file ('phantompeakqualtools/*') from spp_out_mqc.collect()
    file ('phantompeakqualtools/*') from calculateNSCRSC_results.collect()

    output:
    file '*multiqc_report.html'
    file '*multiqc_data'

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    """
    cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
    multiqc -f .
    """
}
