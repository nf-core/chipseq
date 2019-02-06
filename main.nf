#!/usr/bin/env nextflow

/*
========================================================================================
                  C H I P - S E Q   B E S T   P R A C T I C E
========================================================================================
 ChIP-seq Best Practice Analysis Pipeline. Started May 2016.
 #### Homepage / Documentation
 https://github.com/nf-core/chipseq
 @#### Authors
 Chuan Wang <chuan.wang@scilifelab.se>
 Phil Ewels <phil.ewels@scilifelab.se>
 Alex Peltzer <alexander.peltzer@qbic.uni-tuebingen.de>
 Harshil Patel <harshil.patel@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

def helpMessage() {

    log.info nfcoreHeader()
    log.info """

    Usage:

    The typical command for running the pipeline is as follows:

      nextflow run nf-core/chipseq --reads '*_R{1,2}.fastq.gz' --genome GRCh37 --macsconfig 'macssetup.config' -profile docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --fasta                       Path to Fasta reference. Not mandatory when using reference in iGenomes config via --genome
      --macsconfig                  Configuration file for peaking calling using MACS. Format: ChIPSampleID,CtrlSampleID,AnalysisID
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test

    Options:
      --singleEnd                   Specifies that the input is single-end reads

      --allow_multi_align           Secondary alignments and unmapped reads are also reported in addition to primary alignments
      --skipDupRemoval              Skip duplication removal by picard
      --seqCenter                   Text about sequencing center which will be added in the header of output bam files
      --saveAlignedIntermediates    Save the intermediate BAM files from the Alignment step  - not done by default

      --fingerprintBins             Number of genomic bins to use when calculating fingerprint plot. Default: 50000
      --broad                       Run MACS with the --broad flag
      --macsgsize                   Effective genome size for the MACS --gsize option. Should be in the format "2.1e9"
      --saturation                  Run saturation analysis by peak calling with subsets of reads

      --extendReadsLen [int]        Number of base pairs to extend the reads for the deepTools analysis. Default: 100

    References
      --genome                      Name of iGenomes reference
      --bwa_index                   Path to BWA index
      --largeRef                    Build BWA Index for large reference genome (>2Gb)
      --gtf                         Path to GTF file (Ensembl format)
      --bed                         Path to BED file (Ensembl format)
      --blacklist                   Path to blacklist regions (.BED format), used for filtering out called peaks. Note that --blacklist_filtering is required
      --blacklist_filtering         Filter ENCODE blacklisted regions from ChIP-seq peaks. It only works when --genome is set as GRCh37 or GRCm38
      --saveReference               Save the generated reference files in the Results directory

    Trimming options
      --clip_r1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
      --clip_r2 [int]               Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
      --three_prime_clip_r1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed
      --three_prime_clip_r2 [int]   Instructs Trim Galore to re move bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed
      --notrim                      Specifying --notrim will skip the adapter trimming step
      --saveTrimmed                 Save the trimmed Fastq files in the the Results directory

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --rlocation                   Location to save R-libraries used in the pipeline. Default value is ~/R/nxtflow_libs/
      --clusterOptions              Extra SLURM options, used in conjunction with Uppmax.config
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bed = params.genome ? params.genomes[ params.genome ].bed ?: false : false
params.blacklist = params.genome ? params.genomes[ params.genome ].blacklist ?: false : false
params.macsgsize = params.genome ? params.genomes[ params.genome ].macsgsize ?: false : false

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// R library locations
if (params.rlocation){
    nxtflow_libs = file(params.rlocation)
    nxtflow_libs.mkdirs()
}

// Create channels for config files
phantompeakqualtools_mqc_header_ch = Channel.fromPath("$baseDir/assets/phantompeakqualtools_mqc_header", checkIfExists: true)
multiqc_config_ch = Channel.fromPath(params.multiqc_config, checkIfExists: true)
output_docs_ch = Channel.fromPath("$baseDir/docs/output.md", checkIfExists: true)

// Check for file existence
macsconfig = file(params.macsconfig)
if( !macsconfig.exists() ) exit 1, "Missing MACS config: '$macsconfig'. Specify path with --macsconfig"

if( params.bwa_index ){
    bwa_index = Channel
        .fromPath(params.bwa_index, checkIfExists: true)
        .ifEmpty { exit 1, "BWA index not found: ${params.bwa_index}" }
} else if ( params.fasta ){
    fasta = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Fasta file not found: ${params.fasta}" }
} else {
    exit 1, "No reference genome specified!"
}

if( params.gtf ){
    gtf = Channel
        .fromPath(params.gtf, checkIfExists: true)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
}

if( params.bed ){
    bed = Channel
        .fromPath(params.bed, checkIfExists: true)
        .ifEmpty { exit 1, "BED file not found: ${params.bed}" }
}

if( params.blacklist_filtering ){
    blacklist = Channel
        .fromPath(params.blacklist, checkIfExists: true)
        .ifEmpty { exit 1, "Blacklist annotation file not found: ${params.blacklist}" }
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_configvali; raw_reads_fastqc; raw_reads_trimgalore }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_configvali; raw_reads_fastqc; raw_reads_trimgalore }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { raw_reads_configvali; raw_reads_fastqc; raw_reads_trimgalore }
}

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
    .into{ vali_para; macs_para; saturation_para }

// Validate all samples in macs config file
def config_samples = []
for (line in vali_para){
    if (line.getClass().toString() != "class groovyx.gpars.dataflow.operator.PoisonPill") {
        config_samples.add(line[0])
        config_samples.add(line[1])
    }
}
config_samples.removeAll{ it == '' }
config_samples.unique(false)

def fastq_samples = []
for (sample in raw_reads_configvali){
    if (sample.getClass().toString() != "class groovyx.gpars.dataflow.operator.PoisonPill") {
        fastq_samples.add(sample[0].toString() - ~/(.R)?(_R)?(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/)
    }
}

def missing_samples = config_samples - config_samples.intersect(fastq_samples)
if(!missing_samples.isEmpty()){
    exit 1, "No FastQ file found for sample in MACS config: ${missing_samples}"
}

def dropped_samples = fastq_samples - fastq_samples.intersect(config_samples)
if(!dropped_samples.isEmpty()){
    exit 1, "Sample ${dropped_samples} not included in MACS config"
}


log.info nfcoreHeader()
def summary = [:]
summary['Pipeline Name']        = 'nf-core/chipseq'
summary['Pipeline Version']     = workflow.manifest.version
summary['Run Name']             = custom_runName ?: workflow.runName
summary['Reads']                = params.reads
summary['Data Type']            = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']               = params.genome
if(params.bwa_index)  summary['BWA Index'] = params.bwa_index
else if(params.fasta) summary['Fasta File'] = params.fasta
if(params.largeRef)  summary['Build BWA Index for Large Reference'] = params.largeRef
if(params.gtf)  summary['GTF File'] = params.gtf
if(params.bed)  summary['BED File'] = params.bed
summary['Blacklist Filtering']  = params.blacklist_filtering
if(params.blacklist_filtering) summary['Blacklist BED'] = params.blacklist
summary['Save Reference']       = params.saveReference
summary['Multiple Alignments']  = params.allow_multi_align
summary['Duplication Removal']  = params.skipDupRemoval
if(params.seqCenter) summary['Seq Center'] = params.seqCenter
summary['Save Intermeds']       = params.saveAlignedIntermediates
summary['MACS Config']          = params.macsconfig
summary['MACS Broad Peaks']     = params.broad
summary['MACS Genome Size']     = params.macsgsize
summary['Saturation Analysis']  = params.saturation
summary['Extend Reads']         = "$params.extendReadsLen bp"
if( params.notrim ){
    summary['Trimming Step'] = 'Skipped'
} else {
    summary['Trim R1'] = params.clip_r1
    summary['Trim R2'] = params.clip_r2
    summary["Trim 3' R1"] = params.three_prime_clip_r1
    summary["Trim 3' R2"] = params.three_prime_clip_r2
}
summary['Save Trimmed']         = params.saveTrimmed
summary['R libraries']          = params.rlocation
summary['Max Memory']           = params.max_memory
summary['Max CPUs']             = params.max_cpus
summary['Max Time']             = params.max_time
summary['Output Dir']           = params.outdir
summary['Working Dir']          = workflow.workDir
summary['Container Engine']     = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current Home']         = "$HOME"
summary['Current User']         = "$USER"
summary['Current Path']         = "$PWD"
summary['Working Dir']          = workflow.workDir
summary['Output Dir']           = params.outdir
summary['Script Dir']           = workflow.projectDir
summary['Config Profile']       = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']        = params.awsregion
   summary['AWS Queue']         = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "========================================="

// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}

// Check workDir/outdir paths to be S3 buckets if running on AWSBatch
// related: https://github.com/nextflow-io/nextflow/issues/813
if( workflow.profile == 'awsbatch') {
    if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

/*
 * PREPROCESSING - Build BWA index
 */
if(!params.bwa_index && fasta){
    process makeBWAindex {
        tag "$fasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta

        output:
        file "BWAIndex" into bwa_index

        script:
        BWAIndexOption = params.largeRef ? "bwtsw" : 'is'
        """
        bwa index -a $BWAIndexOption $fasta
        mkdir BWAIndex && mv ${fasta}* BWAIndex
        """
    }
}


/*
 * PREPROCESSING - Build BED file
 */
if(!params.bed){
    process makeBED {
        tag "$gtf"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gtf from gtf

        output:
        file "${gtf.baseName}.bed" into bed

        script: // This script is bundled with the pipeline, in nfcore/chipseq/bin/
        """
        gtf2bed $gtf > ${gtf.baseName}.bed
        """
    }
}


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results
    file '.command.out' into fastqc_stdout

    script:
    """
    fastqc -q $reads
    """
}


/*
 * STEP 2 - Trim Galore!
 */
if(params.notrim){
    trimmed_reads = read_files_trimming
    trimgalore_results = []
    trimgalore_fastqc_reports = []
} else {
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
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        if (params.singleEnd) {
            """
            trim_galore --fastqc --gzip $c_r1 $tpc_r1 $reads
            """
        } else {
            """
            trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
            """
        }
    }
}


/*
 * STEP 3.1 - align with bwa
 */
process bwa {
    tag "$prefix"
    publishDir path: { params.saveAlignedIntermediates ? "${params.outdir}/bwa" : params.outdir }, mode: 'copy',
               saveAs: {filename -> params.saveAlignedIntermediates ? filename : null }

    input:
    file reads from trimmed_reads
    file index from bwa_index.collect()

    output:
    file '*.bam' into bwa_bam

    script:
    prefix = reads[0].toString() - ~/(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    filtering = params.allow_multi_align ? '' : "| samtools view -b -q 1 -F 4 -F 256"
    seqCenter = params.seqCenter ? "-R '@RG\\tID:${prefix}\\tCN:${params.seqCenter}'" : ''
    """
    bwa mem -M $seqCenter ${index}/genome.fa $reads | samtools view -bT $index - $filtering > ${prefix}.bam
    """
}


/*
 * STEP 3.2 - post-alignment processing
 */

process samtools {
    tag "${bam.baseName}"
    publishDir path: "${params.outdir}/bwa", mode: 'copy',
               saveAs: { filename ->
                   if (filename.indexOf(".stats.txt") > 0) "stats/$filename"
                   else params.saveAlignedIntermediates ? filename : null
               }

    input:
    file bam from bwa_bam

    output:
    set file("${bam.baseName}.sorted.bam"), file("${bam.baseName}.sorted.bam.bai") into bam_picard, bam_for_mapped
    file "${bam.baseName}.sorted.bed" into bed_total
    file "${bam.baseName}.stats.txt" into samtools_stats

    script:
    """
    samtools sort $bam -o ${bam.baseName}.sorted.bam
    samtools index ${bam.baseName}.sorted.bam
    bedtools bamtobed -i ${bam.baseName}.sorted.bam | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 > ${bam.baseName}.sorted.bed
    samtools stats ${bam.baseName}.sorted.bam > ${bam.baseName}.stats.txt
    """
}


/*
 * STEP 3.3 - Statistics about mapped and unmapped reads against ref genome
 */

process bwa_mapped {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/bwa/mapped", mode: 'copy'

    input:
    set file(bam), file(bai) from bam_for_mapped

    output:
    file 'mapped_refgenome.txt' into bwa_mapped

    script:
    """
    samtools idxstats $bam \\
        | awk -v filename="$bam" '{mapped+=\$3; unmapped+=\$4} END {print filename,"\t",mapped,"\t",unmapped}' \\
        > ${bam.baseName}.mapped_refgenome.txt
    """
}


/*
 * STEP 4 Picard
 */
if (params.skipDupRemoval) {
    bam_picard.into {
        bam_dedup_spp;
        bam_dedup_deepTools;
        bam_dedup_macs;
        bam_dedup_saturation
    }
    picard_reports = Channel.from(false)
} else {
    process picard {
        tag "$prefix"
        publishDir "${params.outdir}/picard", mode: 'copy'

        input:
        set file(bam), file(bai) from bam_picard

        output:
        set file("${prefix}.dedup.sorted.bam"), file("${prefix}.dedup.sorted.bam.bai") into bam_dedup_spp, bam_dedup_deepTools, bam_dedup_macs, bam_dedup_saturation
        file "${prefix}.dedup.sorted.bam" into bed_dedup
        file "${prefix}.picardDupMetrics.txt" into picard_reports

        script:
        prefix = bam[0].toString() - ~/(\.sorted)?(\.bam)?$/
        if( !task.memory ){
            log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
            avail_mem = 3
        } else {
            avail_mem = task.memory.toGiga()
        }
        """
        picard MarkDuplicates \\
            -Xmx${avail_mem}g \\
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
}


/*
 * STEP 5 Read_count_statistics
 */
process countstat {
    tag "${bed[0].baseName}"
    publishDir "${params.outdir}/countstat", mode: 'copy'

    input:
    file bed from params.skipDupRemoval ? bed_total.collect() : bed_total.mix(bed_dedup).collect()

    output:
    file 'read_count_statistics.txt' into countstat_results

    script:
    """
    countstat.pl $bed
    """
}


/*
 * STEP 6.1 Phantompeakqualtools
 */

process phantompeakqualtools {
    tag "$prefix"
    publishDir "${params.outdir}/phantompeakqualtools", mode: 'copy',
                saveAs: {filename -> filename.indexOf(".out") > 0 ? "logs/$filename" : "$filename"}

    input:
    set file(bam), file(bai) from bam_dedup_spp
    file phantompeakqualtools_mqc_header from phantompeakqualtools_mqc_header_ch.collect()

    output:
    file '*.pdf' into spp_plot
    file '*.spp.out' into spp_out, spp_out_mqc
    file '*_mqc.csv' into spp_csv_mqc

    script:
    prefix = bam[0].toString() - ~/(\.dedup)?(\.sorted)?(\.bam)?$/
    """
    run_spp.r -c="$bam" -savp -savd="${prefix}.spp.Rdata" -out="${prefix}.spp.out"
    processSppRdata.r ${prefix}.spp.Rdata ${prefix}.spp.csv
    cat $phantompeakqualtools_mqc_header ${prefix}.spp.csv > ${prefix}_mqc.csv
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
 * STEP 7.1 deepTools bamPEFragmentSize
 */
bam_dedup_deepTools.into {
    bam_dedup_deepTools_bamPEFragmentSize;
    bam_dedup_deepTools_plotFingerprint;
    bam_dedup_deepTools_bamCoverage;
    bam_dedup_deepTools_multiBamSummary
}

process deepTools_bamPEFragmentSize {
    tag "${bam[0].baseName}"
    publishDir "${params.outdir}/deepTools", mode: 'copy'

    input:
    set file(bam), file(bai) from bam_dedup_deepTools_bamPEFragmentSize.collect()

    output:
    file '*.{txt,pdf}' into deepTools_bamPEFragmentSize_results
    file '*.txt' into deepTools_bamPEFragmentSize_multiqc

    when:
    !params.singleEnd

    script:
    """
    bamPEFragmentSize \\
        --binSize 1000 \\
        --bamfiles $bam \\
        --o fragment_size_histogram.pdf \\
        --plotFileFormat pdf \\
        --plotTitle "Paired-end Fragment Size Distribution" \\
        --outRawFragmentLengths bamPEFragmentSize_rawdata.txt
    """
}

/*
 * STEP 7.2 deepTools plotFingerprint
 */
process deepTools_plotFingerprint {
    tag "${bam[0].baseName}"
    publishDir "${params.outdir}/deepTools", mode: 'copy'

    input:
    set file(bam), file(bai) from bam_dedup_deepTools_plotFingerprint.collect()

    output:
    file '*.{txt,pdf}' into deepTools_plotFingerprint_results
    file '*.txt' into deepTools_plotFingerprint_multiqc

    script:
    """
    plotFingerprint \\
        -b $bam \\
        --plotFile fingerprints.pdf \\
        --outRawCounts fingerprint.txt \\
        --extendReads ${params.extendReadsLen} \\
        --skipZeros \\
        --ignoreDuplicates \\
        --numberOfSamples ${params.fingerprintBins} \\
        --binSize 500 \\
        --plotFileFormat pdf \\
        --plotTitle "Fingerprints"
    """
}

/*
 * STEP 7.3 deepTools bamCoverage
 */
process deepTools_bamCoverage {
    tag "${bam[0].baseName}"
    publishDir "${params.outdir}/deepTools", mode: 'copy'

    input:
    set file(bam), file(bai) from bam_dedup_deepTools_bamCoverage.collect()

    output:
    file '*.{bw}' into deepTools_bamCoverage_results
    file '*.txt' into deepTools_bamCoverage_multiqc

    script:
    """
    bamCoverage \\
       -b $bam \\
       --extendReads ${params.extendReadsLen} \\
       --normalizeUsing RPKM \\
       -o ${bam.baseName}.bw
    """
}

/*
 * STEP 7.4 deepTools computeMatrix
 */
process deepTools_computeMatrix {
    tag "${bam[0].baseName}"
    publishDir "${params.outdir}/deepTools", mode: 'copy'

    input:
    file bigwig from deepTools_bamCoverage_results.collect()
    file bed from bed.collect()

    output:
    file '*.{gz}' into deepTools_computeMatrix_results
    file '*.txt' into deepTools_computeMatrix_multiqc

    script:
    """
    computeMatrix \\
        scale-regions \\
        --scoreFileName *.bw \\
        --regionsFileName $bed \\
        --beforeRegionStartLength 3000 \\
        --afterRegionStartLength 3000 \\
        --regionBodyLength 5000 \\
        --outFileName computeMatrix.out.gz \\
        --skipZeros \\
        --smartLabels
    """
}


/*
 * STEP 7.5 deepTools computeMatrix
 */
process deepTools_plotProfile {
    tag "${bam[0].baseName}"
    publishDir "${params.outdir}/deepTools", mode: 'copy'

    input:
    file bigwig from deepTools_computeMatrix_results

    output:
    file '*.{pdf,txt}' into deepTools_plotProfile_results
    file '*.txt' into deepTools_plotProfile_multiqc

    script:
    """
    plotProfile \\
        --matrixFile computeMatrix.out.gz \\
        --outFileName read_distribution_profile.pdf \\
        --plotFileFormat pdf \\
        --outFileNameData read_distribution_profile.txt \\
        --plotTitle "Reads Distribution Profile"
    """
}


/*
 * STEP 7.6 deepTools multiBamSummary
 */
process deepTools_multiBamSummary {
    tag "${bam[0].baseName}"
    publishDir "${params.outdir}/deepTools", mode: 'copy'

    input:
    set file(bam), file(bai) from bam_dedup_deepTools_multiBamSummary.collect()

    output:
    file '*.npz' into deepTools_multiBamSummary_results
    file '*.txt' into deepTools_multiBamSummary_multiqc

    when:
    !(bam instanceof Path)

    script:
    """
    multiBamSummary \\
        bins \\
        --binSize 10000 \\
        --bamfiles $bam \\
        -out multiBamSummary.npz \\
        --extendReads ${params.extendReadsLen} \\
        --ignoreDuplicates \\
        --centerReads \\
        --smartLabels
    """
}


/*
 * STEP 7.7 deepTools plotCorrelation
 */
process deepTools_plotCorrelation {
    tag "${bam[0].baseName}"
    publishDir "${params.outdir}/deepTools", mode: 'copy'

    input:
    file npz from deepTools_multiBamSummary_results.collect()

    output:
    file '*.{pdf,txt}' into deepTools_plotCorrelation_results
    file '*.txt' into deepTools_plotCorrelation_multiqc

    when:
    !(bam instanceof Path)

    script:
    """
    plotCorrelation \\
        -in $npz \\
        -o heatmap_SpearmanCorr.pdf \\
        --plotFileFormat pdf \\
        --outFileCorMatrix heatmap_SpearmanCorr.txt \\
        --corMethod spearman \\
        --skipZeros \\
        --plotTitle "Spearman Correlation of Read Counts" \\
        --whatToPlot heatmap \\
        --colorMap RdYlBu \\
        --plotNumbers
    """
}


/*
 * STEP 7.8 deepTools plotCorrelation
 */
process deepTools_plotPCA {
    tag "${bam[0].baseName}"
    publishDir "${params.outdir}/deepTools", mode: 'copy'

    input:
    file npz from deepTools_multiBamSummary_results.collect()

    output:
    file '*.{pdf,txt}' into deepTools_plotPCA_results
    file '*.txt' into deepTools_plotPCA_multiqc

    when:
    !(bam instanceof Path)

    script:
    """
    plotPCA \\
        -in $npz \\
        -o pcaplot.pdf \\
        --plotFileFormat pdf \\
        --plotTitle "Principal Component Analysis Plot" \\
        --outFileNameData pcaplot.txt \\
        --plotWidth 8
    """
}


/*
 * STEP 8.1 MACS
 */

process macs {
    tag "${analysis_id}"
    publishDir "${params.outdir}/macs", mode: 'copy'

    input:
    set file(bam), file(bai) from bam_dedup_macs.collect()
    set chip_sample_id, ctrl_sample_id, analysis_id from macs_para

    output:
    file '*.{bed,r,narrowPeak}' into macs_results
    file '*.xls' into macs_peaks

    script:
    if(!params.skipDupRemoval){
        chip = "-t ${chip_sample_id}.dedup.sorted.bam"
        ctrl = ctrl_sample_id == '' ? '' : "-c ${ctrl_sample_id}.dedup.sorted.bam"
    }
    else {
        chip = "-t ${chip_sample_id}.sorted.bam"
        ctrl = ctrl_sample_id == '' ? '' : "-c ${ctrl_sample_id}.sorted.bam"
    }
    broad = params.broad ? "--broad" : ''
    """
    macs2 callpeak \\
        $chip \\
        $ctrl \\
        $broad \\
        -f BAM \\
        -g $params.macsgsize \\
        -n $analysis_id \\
        -q 0.01
    """
}


/*
 * STEP 8.2 Saturation analysis
 */
if (params.saturation) {

  process saturation {
     tag "${analysis_id}.${sampling}"
     publishDir "${params.outdir}/macs/saturation", mode: 'copy'

     input:
     set file(bam), file(bai) from bam_dedup_saturation.collect()
     set chip_sample_id, ctrl_sample_id, analysis_id from saturation_para
     each sampling from 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0

     output:
     file '*.xls' into saturation_results

     script:
     if(!params.skipDupRemoval){
         chip_sample = "${chip_sample_id}.dedup.sorted.bam"
         ctrl = ctrl_sample_id == '' ? '' : "-c ${ctrl_sample_id}.dedup.sorted.bam"
     }
     else {
         chip_sample = "${chip_sample_id}.sorted.bam"
         ctrl = ctrl_sample_id == '' ? '' : "-c ${ctrl_sample_id}.sorted.bam"
     }
     broad = params.broad ? "--broad" : ''
     """
     samtools view -b -s ${sampling} ${chip_sample} > ${chip_sample}.${sampling}.bam
     macs2 callpeak \\
         -t ${chip_sample}.${sampling}.bam \\
         $ctrl \\
         $broad \\
         -f BAM \\
         -g $params.macsgsize \\
         -n ${analysis_id}.${sampling} \\
         -q 0.01
     """
  }

  process saturation_r {
     tag "${saturation_results_collection[0].baseName}"
     publishDir "${params.outdir}/macs/saturation", mode: 'copy'

     input:
     file macsconfig from macsconfig
     file countstat from countstat_results
     file saturation_results_collection from saturation_results.collect()

     output:
     file '*.{txt,pdf}' into saturation_summary

     script:
     """
     saturation_results_processing.r $params.rlocation $macsconfig $countstat $saturation_results_collection
     """
  }
}


/*
 * STEP 9 Post peak calling processing
 */

process chippeakanno {
    tag "${macs_peaks_collection[0].baseName}"
    publishDir "${params.outdir}/macs/chippeakanno", mode: 'copy'

    input:
    file macs_peaks_collection from macs_peaks.collect()
    file gtf from gtf

    output:
    file '*.{txt,bed}' into chippeakanno_results

    script:
    filtering = params.blacklist_filtering ? "${params.blacklist}" : "No-filtering"
    """
    post_peak_calling_processing.r $params.rlocation $filtering $gtf $macs_peaks_collection
    """
}

/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    trim_galore --version > v_trim_galore.txt
    echo \$(bwa 2>&1) > v_bwa.txt
    samtools --version > v_samtools.txt
    bedtools --version > v_bedtools.txt
    picard MarkDuplicates --version &> v_picard.txt  || true
    echo \$(plotFingerprint --version 2>&1) > v_deeptools.txt
    echo \$(macs2 --version 2>&1) > v_macs2.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


/*
 * STEP 10 MultiQC
 */

process multiqc {
    tag "$prefix"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config from multiqc_config_ch.collect()
    file (fastqc:'fastqc/*') from fastqc_results.collect()
    file ('trimgalore/*') from trimgalore_results.collect()
    file ('samtools/*') from samtools_stats.collect()
    file ('picard/*') from picard_reports.collect()
    file ('deeptools/bamPEFragmentSize/*') from deepTools_bamPEFragmentSize_multiqc.collect()
    file ('deeptools/plotFingerprint/*') from deepTools_plotFingerprint_multiqc.collect()
    file ('deeptools/bamCoverage/*') from deepTools_bamCoverage_multiqc.collect()
    file ('deeptools/computeMatrix/*') from deepTools_computeMatrix_multiqc.collect()
    file ('deeptools/plotProfile/*') from deepTools_plotProfile_multiqc.collect()
    file ('deeptools/multiBamSummary/*') from deepTools_multiBamSummary_multiqc.collect()
    file ('deeptools/plotCorrelation/*') from deepTools_plotCorrelation_multiqc.collect()
    file ('deeptools/plotPCA/*') from deepTools_plotPCA_multiqc.collect()
    file ('phantompeakqualtools/*') from spp_out_mqc.collect()
    file ('phantompeakqualtools/*') from spp_csv_mqc.collect()
    file ('software_versions/*') from software_versions_yaml.collect()

    output:
    file '*multiqc_report.html' into multiqc_report
    file '*_data' into multiqc_data
    file '.command.err' into multiqc_stderr

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config \\
        -m custom_content -m fastqc -m cutadapt -m samtools -m picard -m deeptools -m phantompeakqualtools
    """
}

/*
 * STEP 11 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs from output_docs_ch

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/chipseq] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/chipseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/chipseq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/chipseq] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/chipseq] Pipeline Complete"

}

def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """
    ${c_dim}====================================================${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/chipseq v${workflow.manifest.version}${c_reset}
    ${c_dim}====================================================${c_reset}
    """.stripIndent()
}
