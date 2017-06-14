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
*/



/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 1.2

// Configurable variables
params.name = false
params.project = false
params.genome = false
params.genomes = []
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
params.reads = "data/*{1,2}*.fastq.gz"
params.macsconfig = "data/macsconfig"
params.extendReadsLen = 100
params.notrim = false
params.saveReference = false
params.saveTrimmed = false
params.saveAlignedIntermediates = false
params.broad = false
params.outdir = './results'
params.email = false

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

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Custom trimming options
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

// Config / R locations
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
multiqc_config = file(params.multiqc_config)
params.rlocation = false
if (params.rlocation){
    nxtflow_libs = file(params.rlocation)
    nxtflow_libs.mkdirs()
}


/*
 * Create a channel for input read files
 */
params.singleEnd = false
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nIf this is single-end data, please specify --singleEnd on the command line." }
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


// Header log info
log.info "========================================="
log.info " NGI-ChIPseq: ChIP-Seq Best Practice v${version}"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']       = params.genome
if(params.bwa_index)  summary['BWA Index'] = params.bwa_index
else if(params.fasta) summary['Fasta Ref'] = params.fasta
summary['MACS Config']    = params.macsconfig
summary['MACS broad peaks'] = params.broad
summary['Extend Reads']   = "$params.extendReadsLen bp"
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['R libraries']    = params.rlocation
summary['Script dir']     = workflow.projectDir
summary['Save Reference'] = params.saveReference
summary['Save Trimmed']   = params.saveTrimmed
summary['Save Intermeds'] = params.saveAlignedIntermediates
if(params.notrim)       summary['Trimming Step'] = 'Skipped'
if( params.clip_r1 > 0) summary['Trim R1'] = params.clip_r1
if( params.clip_r2 > 0) summary['Trim R2'] = params.clip_r2
if( params.three_prime_clip_r1 > 0) summary["Trim 3' R1"] = params.three_prime_clip_r1
if( params.three_prime_clip_r2 > 0) summary["Trim 3' R2"] = params.three_prime_clip_r2
summary['Config Profile'] = (workflow.profile == 'standard' ? 'UPPMAX' : workflow.profile)
if(params.project) summary['UPPMAX Project'] = params.project
if(params.email) summary['E-mail Address'] = params.email
if(workflow.commitId) summary['Pipeline Commit']= workflow.commitId
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "===================================="


/*
 * PREPROCESSING - Build BWA index
 */
if(!params.bwa_index && fasta){
    process makeBWAindex {
        tag fasta
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

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
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

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
    publishDir path: { params.saveAlignedIntermediates ? "${params.outdir}/bwa" : params.outdir }, mode: 'copy',
               saveAs: {filename -> params.saveAlignedIntermediates ? filename : null }

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
        log.warn "[Picard MarkDuplicates] Available memory not known - defaulting to 6GB ($prefix)"
        avail_mem = 6000
    } else {
        avail_mem = task.memory.toMega()
        if( avail_mem <= 0){
            avail_mem = 6000
            log.warn "[Picard MarkDuplicates] Available memory 0 - defaulting to 6GB ($prefix)"
        } else if( avail_mem < 250){
            avail_mem = 250
            log.warn "[Picard MarkDuplicates] Available memory under 250MB - defaulting to 250MB ($prefix)"
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
    publishDir "${params.outdir}/phantompeakqualtools", mode: 'copy',
                saveAs: {filename -> filename.indexOf(".out") > 0 ? "logs/$filename" : "$filename"}

    input:
    file bam from bam_dedup_spp

    output:
    file '*.pdf' into spp_results
    file '*.spp.out' into spp_out, spp_out_mqc

    script:
    prefix = bam[0].toString() - ~/(\.dedup)?(\.sorted)?(\.bam)?$/
    """
    SCRIPT_PATH=\$(which run_spp.R)
    Rscript \$SCRIPT_PATH -c="$bam" -savp -out="${prefix}.spp.out"
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
    broad = params.broad ? "--broad" : ''
    """
    macs2 callpeak \\
        -t ${chip_sample_id}.dedup.sorted.bam \\
        $ctrl \\
        $broad \\
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
    file multiqc_config
    file (fastqc:'fastqc/*') from fastqc_results.collect()
    file ('trimgalore/*') from trimgalore_results.collect()
    file ('bwa/*') from bwa_logs.collect()
    file ('picard/*') from picard_reports.collect()
    file ('phantompeakqualtools/*') from spp_out_mqc.collect()
    file ('phantompeakqualtools/*') from calculateNSCRSC_results.collect()

    output:
    file '*multiqc_report.html' into multiqc_report
    file '*_data'
    val prefix into multiqc_prefix

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config . 2>&1
    """
}

/*
 * STEP 11 - Output Description HTML
 */
process output_documentation {
    tag "$prefix"
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    val prefix from multiqc_prefix

    output:
    file "results_description.html"

    script:
    def rlocation = params.rlocation ?: ''
    """
    markdown_to_html.r $baseDir/docs/output.md results_description.html $rlocation
    """
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Build the e-mail subject and header
    def subject = "NGI-ChIPseq Pipeline Complete: $workflow.runName"
    subject += "\nContent-Type: text/html"

    // Set up the e-mail variables
    def email_fields = [:]
    email_fields['version'] = version
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
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container

    // Render the e-mail HTML template
    def f = new File("$baseDir/assets/summary_email.html")
    def engine = new groovy.text.GStringTemplateEngine()
    def template = engine.createTemplate(f).make(email_fields)
    def email_html = template.toString()

    // Send the HTML e-mail
    if (params.email) {
        [ 'mail', '-s', subject, params.email ].execute() << email_html
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_f = new File( output_d, "pipeline_report.html" )
    output_f.withWriter { w ->
        w << email_html
    }

}
