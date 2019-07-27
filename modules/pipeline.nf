///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                  SET DEFAULT PARAMETER VALUES                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Options: Generic
params.design = false
params.singleEnd = false
params.seq_center = false
params.fragment_size = 200
params.fingerprint_bins = 500000

// Options: References
params.genome = false
params.fasta = false
params.gtf = false
params.gene_bed = false
params.tss_bed = false
params.bwa_index = false
params.blacklist = false
params.macs_gsize = false
params.saveGenomeIndex = false

// Options: Trimming
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0
params.skipTrimming = false
params.saveTrimmed = false

// Options: Alignments
params.keepDups = false
params.keepMultiMap = false
params.saveAlignedIntermediates = false

// Options: Peaks
params.narrowPeak = false
params.broad_cutoff = 0.1
params.min_reps_consensus = 1
params.saveMACSPileup = false
params.skipDiffAnalysis = false

// Options: QC
params.skipFastQC = false
params.skipPicardMetrics = false
params.skipPreseq = false
params.skipPlotProfile = false
params.skipPlotFingerprint = false
params.skipSpp = false
params.skipIGV = false
params.skipMultiQC = false
params.skipMultiQCStats = false

// Options: AWSBatch
params.awsqueue = false
params.awsregion = 'eu-west-1'

// Options: Custom config
params.config_profile_description = false
params.config_profile_contact = false
params.config_profile_url = false

// Options: Other
params.outdir = './results'
params.run_name = false
params.email = false
params.monochrome_logs = false
params.hostnames = false
params.maxMultiqcEmailFileSize = 25.MB
params.max_memory = 128.GB
params.max_cpus = 16
params.max_time = 240.h

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         PRINT HELP MESSAGE                          -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

def print_help() {
    log.info nfcore_header()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

      nextflow run nf-core/chipseq --design design.csv --genome GRCh37 -profile docker

    Mandatory arguments:
      --design                      Comma-separated file containing information about the samples in the experiment (see docs/usage.md)
      --fasta                       Path to Fasta reference. Not mandatory when using reference in iGenomes config via --genome
      --gtf                         Path to GTF file in Ensembl format. Not mandatory when using reference in iGenomes config via --genome
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test

    Generic
      --singleEnd                   Specifies that the input is single-end reads
      --seq_center                  Sequencing center information to be added to read group of BAM files
      --fragment_size [int]         Estimated fragment size used to extend single-end reads (Default: 200)
      --fingerprint_bins [int]      Number of genomic bins to use when calculating fingerprint plot (Default: 500000)

    References                      If not specified in the configuration file or you wish to overwrite any of the references
      --genome                      Name of iGenomes reference
      --bwa_index                   Full path to directory containing BWA index including base name i.e. /path/to/index/genome.fa
      --gene_bed                    Path to BED file containing gene intervals
      --tss_bed                     Path to BED file containing transcription start sites
      --macs_gsize                  Effective genome size parameter required by MACS2. If using iGenomes config, values have only been provided when --genome is set as GRCh37, GRCm38, hg19, mm10, BDGP6 and WBcel235
      --blacklist                   Path to blacklist regions (.BED format), used for filtering alignments
      --saveGenomeIndex             If generated by the pipeline save the BWA index in the results directory

    Trimming
      --clip_r1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads) (Default: 0)
      --clip_r2 [int]               Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only) (Default: 0)
      --three_prime_clip_r1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed (Default: 0)
      --three_prime_clip_r2 [int]   Instructs Trim Galore to re move bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed (Default: 0)
      --skipTrimming                Skip the adapter trimming step
      --saveTrimmed                 Save the trimmed FastQ files in the the results directory

    Alignments
      --keepDups                    Duplicate reads are not filtered from alignments
      --keepMultiMap                Reads mapping to multiple locations are not filtered from alignments
      --saveAlignedIntermediates    Save the intermediate BAM files from the alignment step - not done by default

    Peaks
      --narrowPeak                  Run MACS2 in narrowPeak mode
      --broad_cutoff [float]        Specifies broad cutoff value for MACS2. Only used when --narrowPeak isnt specified (Default: 0.1)
      --min_reps_consensus          Number of biological replicates required from a given condition for a peak to contribute to a consensus peak (Default: 1)
      --saveMACSPileup              Instruct MACS2 to create bedGraph files normalised to signal per million reads
      --skipDiffAnalysis            Skip differential binding analysis

    QC
      --skipFastQC                  Skip FastQC
      --skipPicardMetrics           Skip Picard CollectMultipleMetrics
      --skipPreseq                  Skip Preseq
      --skipPlotProfile             Skip deepTools plotProfile
      --skipPlotFingerprint         Skip deepTools plotFingerprint
      --skipSpp                     Skip Phantompeakqualtools
      --skipIGV                     Skip IGV
      --skipMultiQC                 Skip MultiQC
      --skipMultiQCStats            Exclude general statistics table from MultiQC report

    Other
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --maxMultiqcEmailFileSize     Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                STORE AND PRINT PARAMETER VALUES                     -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

def create_summary() {
    log.info nfcore_header()
    def summary = [:]
    summary['Run Name']             = params.run_name ?: workflow.runName
    summary['Data Type']            = params.singleEnd ? 'Single-End' : 'Paired-End'
    summary['Design File']          = params.design
    summary['Genome']               = params.genome ?: 'Not supplied'
    summary['Fasta File']           = params.fasta
    summary['GTF File']             = params.gtf
    if (params.gene_bed)            summary['Gene BED File'] = params.gene_bed
    if (params.tss_bed)             summary['TSS BED File'] = params.tss_bed
    if (params.bwa_index)           summary['BWA Index'] = params.bwa_index
    if (params.blacklist)           summary['Blacklist BED'] = params.blacklist
    summary['MACS2 Genome Size']    = params.macs_gsize ?: 'Not supplied'
    summary['Min Consensus Reps']   = params.min_reps_consensus
    if (params.macs_gsize)          summary['MACS2 Narrow Peaks'] = params.narrowPeak ? 'Yes' : 'No'
    if (!params.narrowPeak)         summary['MACS2 Broad Cutoff'] = params.broad_cutoff
    if (params.skipTrimming){
        summary['Trimming Step']    = 'Skipped'
    } else {
        summary['Trim R1']          = "$params.clip_r1 bp"
        summary['Trim R2']          = "$params.clip_r2 bp"
        summary["Trim 3' R1"]       = "$params.three_prime_clip_r1 bp"
        summary["Trim 3' R2"]       = "$params.three_prime_clip_r2 bp"
    }
    if (params.seq_center)          summary['Sequencing Center'] = params.seq_center
    if (params.singleEnd)           summary['Fragment Size'] = "$params.fragment_size bp"
    summary['Fingerprint Bins']     = params.fingerprint_bins
    if (params.keepDups)            summary['Keep Duplicates'] = 'Yes'
    if (params.keepMultiMap)        summary['Keep Multi-mapped'] = 'Yes'
    summary['Save Genome Index']    = params.saveGenomeIndex ? 'Yes' : 'No'
    if (params.saveTrimmed)         summary['Save Trimmed'] = 'Yes'
    if (params.saveAlignedIntermediates) summary['Save Intermeds'] =  'Yes'
    if (params.saveMACSPileup)      summary['Save MACS2 Pileup'] = 'Yes'
    if (params.skipDiffAnalysis)    summary['Skip Diff Analysis'] = 'Yes'
    if (params.skipFastQC)          summary['Skip FastQC'] = 'Yes'
    if (params.skipPicardMetrics)   summary['Skip Picard Metrics'] = 'Yes'
    if (params.skipPreseq)          summary['Skip Preseq'] = 'Yes'
    if (params.skipPlotProfile)     summary['Skip plotProfile'] = 'Yes'
    if (params.skipPlotFingerprint) summary['Skip plotFingerprint'] = 'Yes'
    if (params.skipSpp)             summary['Skip spp'] = 'Yes'
    if (params.skipIGV)             summary['Skip IGV'] = 'Yes'
    if (params.skipMultiQC)         summary['Skip MultiQC'] = 'Yes'
    if (params.skipMultiQCStats)    summary['Skip MultiQC Stats'] = 'Yes'
    summary['Max Resources']        = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
    if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
    summary['Output Dir']           = params.outdir
    summary['Launch Dir']           = workflow.launchDir
    summary['Working Dir']          = workflow.workDir
    summary['Script Dir']           = workflow.projectDir
    summary['User']                 = workflow.userName
    if (workflow.profile == 'awsbatch'){
       summary['AWS Region']        = params.awsregion
       summary['AWS Queue']         = params.awsqueue
    }
    summary['Config Profile']       = workflow.profile
    if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
    if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
    if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
    if(params.email) {
      summary['E-mail Address']     = params.email
      summary['MultiQC Max Size']   = params.maxMultiqcEmailFileSize
    }
    log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
    log.info "\033[2m----------------------------------------------------\033[0m"

    return summary
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                      PIPELINE PROCESSES                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * Show a big warning message if we're not running MACS
 */
if (!params.macs_gsize){
    def warnstring = params.genome ? "supported for '${params.genome}'" : 'supplied'
    log.warn "=================================================================\n" +
             "  WARNING! MACS genome size parameter not $warnstring.\n" +
             "  Peak calling, annotation and differential analysis will be skipped.\n" +
             "  Please specify value for '--macs_gsize' to run these steps.\n" +
             "======================================================================="
}

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".csv") > 0) filename
            else null
        }

    output:
    file 'software_versions_mqc.yaml'
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    trim_galore --version > v_trim_galore.txt
    echo \$(bwa 2>&1) > v_bwa.txt
    samtools --version > v_samtools.txt
    bedtools --version > v_bedtools.txt
    echo \$(bamtools --version 2>&1) > v_bamtools.txt
    echo \$(plotFingerprint --version 2>&1) > v_deeptools.txt || true
    picard MarkDuplicates --version &> v_picard.txt  || true
    echo \$(R --version 2>&1) > v_R.txt
    python -c "import pysam; print(pysam.__version__)" > v_pysam.txt
    echo \$(macs2 --version 2>&1) > v_macs2.txt
    touch v_homer.txt
    echo \$(featureCounts -v 2>&1) > v_featurecounts.txt
    preseq &> v_preseq.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


/*
 * Output description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}


/*
 * Create summary file for MultiQC
 */
def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-chipseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/chipseq Workflow Summary'
    section_href: 'https://github.com/nf-core/chipseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         COMPLETION EMAIL                            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * Completion e-mail notification
 */
def send_email(summary) {

    // Set up the e-mail variables
    def subject = "[nf-core/chipseq] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/chipseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = params.run_name ?: workflow.runName
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
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // // // On success try attach the multiqc report
    // // def mqc_report = null
    // // try {
    // //     if (workflow.success) {
    // //         mqc_report = ch_multiqc_report.getVal()
    // //         if (mqc_report.getClass() == ArrayList){
    // //             log.warn "[nf-core/chipseq] Found multiple reports from process 'multiqc', will use only one"
    // //             mqc_report = mqc_report[0]
    // //         }
    // //     }
    // // } catch (all) {
    // //     log.warn "[nf-core/chipseq] Could not attach MultiQC report to summary email"
    // // }
    //
    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // // Render the sendmail template
    // //def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: false, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes()]
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
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
        output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if(workflow.success){
        log.info "${c_purple}[nf-core/chipseq]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        check_host_name()
        log.info "${c_purple}[nf-core/chipseq]${c_red} Pipeline completed with errors${c_reset}"
    }

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       TEMPLATE CODE                                 -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * AWS batch settings
 */
if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

/*
 * Print nf-core header
 */
def nfcore_header(){
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

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/chipseq v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}


/*
 * Check host name if provided on nf-core/configs
 */
def check_host_name(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if(params.hostnames){
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if(hostname.contains(hname) && !workflow.profile.contains(prof)){
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}