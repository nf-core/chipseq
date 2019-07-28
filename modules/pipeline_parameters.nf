
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
params.maxMultiqcEmailFileSize = 25.MB
params.max_memory = 128.GB
params.max_cpus = 16
params.max_time = 240.h

include 'nfcore_header' params(params)


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
 * Print help
 */
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


/*
 * Print parameter summary
 */
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
