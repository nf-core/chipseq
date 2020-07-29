#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/chipseq
========================================================================================
 nf-core/chipseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/chipseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
 * Print help message if required
 */
if (params.help) {
    def command = "nextflow run nf-core/chipseq --input design.csv --genome GRCh37 -profile docker"
    log.info Headers.nf_core(workflow, params.monochrome_logs)
    log.info Schema.params_help("$baseDir/nextflow_schema.json", command)
    exit 0
}

////////////////////////////////////////////////////
/* --         DEFAULT PARAMETER VALUES         -- */
////////////////////////////////////////////////////

/*
 * Reference genomes
 */
// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Configurable variables
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.gene_bed = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.macs_gsize = params.genome ? params.genomes[ params.genome ].macs_gsize ?: false : false
params.blacklist = params.genome ? params.genomes[ params.genome ].blacklist ?: false : false
params.anno_readme = params.genome ? params.genomes[ params.genome ].readme ?: false : false

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

/*
 * Validate parameters
 */
if (params.input)     { ch_input = file(params.input, checkIfExists: true) } else { exit 1, 'Samples design file not specified!' }
if (params.gtf)       { ch_gtf = file(params.gtf, checkIfExists: true) } else { exit 1, 'GTF annotation file not specified!' }
if (params.gene_bed)  { ch_gene_bed = file(params.gene_bed, checkIfExists: true) }
if (params.blacklist) { ch_blacklist = Channel.fromPath(params.blacklist, checkIfExists: true) } else { ch_blacklist = Channel.empty() }

if (params.fasta) {
    ch_fasta = file(params.fasta, checkIfExists: true)
} else {
    exit 1, 'Fasta file not specified!'
}

// Save AWS IGenomes file containing annotation version
if (params.anno_readme && file(params.anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(params.anno_readme).copyTo("${params.outdir}/genome/")
}

// If --gtf is supplied along with --genome
// Make gene bed from supplied --gtf instead of using iGenomes one automatically
def makeBED = false
if (!params.gene_bed) {
    makeBED = true
} else if (params.genome && params.gtf) {
    if (params.genomes[ params.genome ].gtf != params.gtf) {
        makeBED = true
    }
}

/*
 * Check parameters
 */
Checks.aws_batch(workflow, params)     // Check AWS batch settings
Checks.hostname(workflow, params, log) // Check the hostnames against configured profiles
Checks.macs2_warn(params, log)         // Show a big warning message if we're not running MACS

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

/*
 * Stage config files
 */
ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)

// JSON files required by BAMTools for alignment filtering
ch_bamtools_filter_se_config = file(params.bamtools_filter_se_config, checkIfExists: true)
ch_bamtools_filter_pe_config = file(params.bamtools_filter_pe_config, checkIfExists: true)
//ch_bamtools_filter_config

// Header files for MultiQC
ch_peak_count_header = file("$baseDir/assets/multiqc/peak_count_header.txt", checkIfExists: true)
ch_frip_score_header = file("$baseDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)
ch_peak_annotation_header = file("$baseDir/assets/multiqc/peak_annotation_header.txt", checkIfExists: true)
ch_deseq2_pca_header = file("$baseDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_deseq2_clustering_header = file("$baseDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
ch_spp_correlation_header = file("$baseDir/assets/multiqc/spp_correlation_header.txt", checkIfExists: true)
ch_spp_nsc_header = file("$baseDir/assets/multiqc/spp_nsc_header.txt", checkIfExists: true)
ch_spp_rsc_header = file("$baseDir/assets/multiqc/spp_rsc_header.txt", checkIfExists: true)

////////////////////////////////////////////////////
/* --          PARAMETER SUMMARY               -- */
////////////////////////////////////////////////////

/*
 * Print parameter summary
 */
// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
run_name = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    run_name = workflow.runName
}
summary = Schema.params_summary(workflow, params, run_name)
log.info Headers.nf_core(workflow, params.monochrome_logs)
log.info summary.collect { k,v -> "${k.padRight(22)}: $v" }.join("\n")
log.info "-\033[2m----------------------------------------------------\033[0m-"

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { GTF2BED                  } from './modules/local/process/gtf2bed'
include { GET_CHROM_SIZES          } from './modules/local/process/get_chrom_sizes'
include { MAKE_GENOME_FILTER       } from './modules/local/process/make_genome_filter'
include { BEDTOOLS_GENOMECOV       } from './modules/local/process/bedtools_genomecov'
include { PLOT_HOMER_ANNOTATEPEAKS } from './modules/local/process/plot_homer_annotatepeaks'
include { PLOT_MACS2_QC            } from './modules/local/process/plot_macs2_qc'
include { MACS2_CONSENSUS          } from './modules/local/process/macs2_consensus'
// include { FRIP_SCORE               } from './modules/local/process/frip_score'
include { OUTPUT_DOCUMENTATION     } from './modules/local/process/output_documentation'
include { GET_SOFTWARE_VERSIONS    } from './modules/local/process/get_software_versions'
include { MULTIQC                  } from './modules/local/process/multiqc'

include { INPUT_CHECK              } from './modules/local/subworkflow/input_check'
include { BAM_CLEAN                } from './modules/local/subworkflow/bam_clean'

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

include { BWA_INDEX                     } from './modules/nf-core/software/bwa_index'
include { PICARD_MERGESAMFILES          } from './modules/nf-core/software/picard_mergesamfiles'
include { PICARD_COLLECTMULTIPLEMETRICS } from './modules/nf-core/software/picard_collectmultiplemetrics'
include { PRESEQ_LCEXTRAP               } from './modules/nf-core/software/preseq_lcextrap'
include { UCSC_BEDRAPHTOBIGWIG          } from './modules/nf-core/software/ucsc_bedgraphtobigwig'
include { DEEPTOOLS_COMPUTEMATRIX       } from './modules/nf-core/software/deeptools_computematrix'
include { DEEPTOOLS_PLOTPROFILE         } from './modules/nf-core/software/deeptools_plotprofile'
include { DEEPTOOLS_PLOTHEATMAP         } from './modules/nf-core/software/deeptools_plotheatmap'
include { DEEPTOOLS_PLOTFINGERPRINT     } from './modules/nf-core/software/deeptools_plotfingerprint'
include { PHANTOMPEAKQUALTOOLS          } from './modules/nf-core/software/phantompeakqualtools'
include { MACS2_CALLPEAK                } from './modules/nf-core/software/macs2_callpeak'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_MACS2
          HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_CONSENSUS } from './modules/nf-core/software/homer_annotatepeaks'
//include { SUBREAD_FEATURECOUNTS       } from './modules/nf-core/software/subread_featurecounts'

include { FASTQC_TRIMGALORE             } from './modules/nf-core/subworkflow/fastqc_trimgalore'
include { MAP_BWA_MEM                   } from './modules/nf-core/subworkflow/map_bwa_mem'
include { MARK_DUPLICATES_PICARD        } from './modules/nf-core/subworkflow/mark_duplicates_picard'

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {

    /*
     * Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK (
        ch_input,
        params.seq_center,
        params.modules['samplesheet_check']
    )

    /*
     * Prepare genome files
     */
    ch_index = params.bwa_index ? Channel.value(file(params.bwa_index)) : BWA_INDEX ( ch_fasta, params.modules['bwa_index'] ).index

    if (makeBED) { ch_gene_bed = GTF2BED ( ch_gtf, params.modules['gtf2bed'] ) }

    MAKE_GENOME_FILTER (
        GET_CHROM_SIZES ( ch_fasta, params.modules['get_chrom_sizes'] ).sizes,
        ch_blacklist.ifEmpty([]),
        params.modules['make_genome_filter']
    )
    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(MAKE_GENOME_FILTER.out.version.first().ifEmpty(null))

    /*
     * Read QC & trimming
     */
    nextseq = params.trim_nextseq > 0 ? " --nextseq ${params.trim_nextseq}" : ''
    params.modules['trimgalore'].args += nextseq
    FASTQC_TRIMGALORE (
        INPUT_CHECK.out.reads,
        params.skip_fastqc,
        params.skip_trimming,
        params.modules['fastqc'],
        params.modules['trimgalore']
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.fastqc_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.trimgalore_version.first().ifEmpty(null))

    /*
     * Map reads & BAM QC
     */
    score = params.bwa_min_score ? " -T ${params.bwa_min_score}" : ''
    params.modules['bwa_mem'].args += score
    MAP_BWA_MEM (
        FASTQC_TRIMGALORE.out.reads,
        ch_index,
        ch_fasta,
        params.modules['bwa_mem'],
        params.modules['samtools_sort_lib']
    )
    ch_software_versions = ch_software_versions.mix(MAP_BWA_MEM.out.bwa_version.first())
    ch_software_versions = ch_software_versions.mix(MAP_BWA_MEM.out.samtools_version.first().ifEmpty(null))

    /*
     * Merge resequenced BAM files
     */
    MAP_BWA_MEM
        .out
        .bam
        .map {
            meta, bam ->
                fmeta = meta.findAll { it.key != 'read_group' }
                fmeta.id = fmeta.id.split('_')[0..-2].join('_')
                [ fmeta, bam ] }
       .groupTuple(by: [0])
       .map { it ->  [ it[0], it[1].flatten() ] }
       .set { ch_sort_bam }

    PICARD_MERGESAMFILES (
        ch_sort_bam,
        params.modules['picard_mergesamfiles']
    )
    ch_software_versions = ch_software_versions.mix(PICARD_MERGESAMFILES.out.version.first().ifEmpty(null))

    /*
     * Mark duplicates & filter BAM files
     */
    MARK_DUPLICATES_PICARD (
        PICARD_MERGESAMFILES.out.bam,
        params.modules['picard_markduplicates'],
        params.modules['samtools_sort_merged_lib']
    )

    // Fix getting name sorted BAM here for PE/SE
    BAM_CLEAN (
        MARK_DUPLICATES_PICARD.out.bam.join(MARK_DUPLICATES_PICARD.out.bai, by: [0]),
        MAKE_GENOME_FILTER.out.bed.collect(),
        ch_bamtools_filter_se_config,
        ch_bamtools_filter_pe_config,
        params.modules['bam_filter'],
        params.modules['bam_remove_orphans'],
        params.modules['samtools_sort_filter']
    )
    ch_software_versions = ch_software_versions.mix(BAM_CLEAN.out.bamtools_version.first().ifEmpty(null))

    /*
     * Post alignment QC
     */
    PICARD_COLLECTMULTIPLEMETRICS (
        BAM_CLEAN.out.bam,
        ch_fasta,
        params.modules['picard_collectmultiplemetrics']
    )

    PRESEQ_LCEXTRAP (
        BAM_CLEAN.out.bam,
        params.modules['preseq_lcextrap']
    )
    ch_software_versions = ch_software_versions.mix(PRESEQ_LCEXTRAP.out.version.first().ifEmpty(null))

    PHANTOMPEAKQUALTOOLS (
        BAM_CLEAN.out.bam,
        params.modules['phantompeakqualtools']
    )
    ch_software_versions = ch_software_versions.mix(PHANTOMPEAKQUALTOOLS.out.version.first().ifEmpty(null))

    /*
     * Coverage tracks
     */
    BEDTOOLS_GENOMECOV (
        BAM_CLEAN.out.bam.join(BAM_CLEAN.out.flagstat, by: [0]),
        params.modules['bedtools_genomecov']
    )

    UCSC_BEDRAPHTOBIGWIG (
        BEDTOOLS_GENOMECOV.out.bedgraph,
        GET_CHROM_SIZES.out.sizes,
        params.modules['ucsc_bedgraphtobigwig']
    )
    ch_software_versions = ch_software_versions.mix(UCSC_BEDRAPHTOBIGWIG.out.version.first().ifEmpty(null))

    /*
     * Coverage plots
     */
    DEEPTOOLS_COMPUTEMATRIX (
        UCSC_BEDRAPHTOBIGWIG.out.bigwig,
        ch_gene_bed,
        params.modules['deeptools_computematrix']
    )
    ch_software_versions = ch_software_versions.mix(DEEPTOOLS_COMPUTEMATRIX.out.version.first().ifEmpty(null))

    DEEPTOOLS_PLOTPROFILE (
        DEEPTOOLS_COMPUTEMATRIX.out.matrix,
        params.modules['deeptools_plotprofile']
    )

    DEEPTOOLS_PLOTHEATMAP (
        DEEPTOOLS_COMPUTEMATRIX.out.matrix,
        params.modules['deeptools_plotheatmap']
    )

    /*
     * Refactor channels: [ val(meta), [ ip_bam, control_bam ] [ ip_bai, control_bai ] ]
     */
    BAM_CLEAN
        .out
        .bam
        .join ( BAM_CLEAN.out.bai, by: [0] )
        .map { meta, bam, bai -> meta.control ? null : [ meta.id, [ bam ] , [ bai ] ] }
        .set { ch_ip_control_bam }

    BAM_CLEAN
        .out
        .bam
        .join ( BAM_CLEAN.out.bai, by: [0] )
        .map { meta, bam, bai -> meta.control ? [ meta.control, meta, [ bam ], [ bai ] ] : null }
        .combine(ch_ip_control_bam, by: 0)
        .map { it -> [ it[1] , it[2] + it[4], it[3] + it[5] ] }
        .set { ch_ip_control_bam }

    /*
     * plotFingerprint for IP and control together
     */
    params.modules['deeptools_plotfingerprint'].args += " --numberOfSamples $params.fingerprint_bins"
    DEEPTOOLS_PLOTFINGERPRINT (
        ch_ip_control_bam,
        params.modules['deeptools_plotfingerprint']
    )

    if (params.macs_gsize) {

        /*
         * Call peaks
         */
        peakType = params.narrow_peak ? 'narrowPeak' : 'broadPeak'
        broad = params.narrow_peak ? '' : "--broad --broad-cutoff ${params.broad_cutoff}"
        pileup = params.save_macs_pileup ? '--bdg --SPMR' : ''
        fdr = params.macs_fdr ? "--qvalue ${params.macs_fdr}" : ''
        pvalue = params.macs_pvalue ? "--pvalue ${params.macs_pvalue}" : ''
        params.modules['macs2_callpeak'].publish_dir += "/$peakType"
        params.modules['macs2_callpeak'].args += " $broad $pileup $fdr $pvalue"

        // Create channel: [ val(meta), ip_bam, control_bam ]
        ch_ip_control_bam
            .map { it -> [ it[0] , it[1][0], it[1][1] ] }
            .set { ch_ip_control_bam }

        MACS2_CALLPEAK (
            ch_ip_control_bam,
            params.macs_gsize,
            params.modules['macs2_callpeak']
        )
        ch_software_versions = ch_software_versions.mix(MACS2_CALLPEAK.out.version.first().ifEmpty(null))

        params.modules['plot_macs2_qc'].publish_dir += "/$peakType/qc"
        PLOT_MACS2_QC (
            MACS2_CALLPEAK.out.peak.collect{it[1]},
            params.modules['plot_macs2_qc']
        )

        params.modules['homer_annotatepeaks_macs2'].publish_dir += "/$peakType"
        HOMER_ANNOTATEPEAKS_MACS2 (
            MACS2_CALLPEAK.out.peak,
            ch_fasta,
            ch_gtf,
            params.modules['homer_annotatepeaks_macs2']
        )
        ch_software_versions = ch_software_versions.mix(HOMER_ANNOTATEPEAKS_MACS2.out.version.first().ifEmpty(null))

        params.modules['plot_homer_annotatepeaks'].publish_dir += "/$peakType/qc"
        PLOT_HOMER_ANNOTATEPEAKS (
            HOMER_ANNOTATEPEAKS_MACS2.out.txt.collect{it[1]},
            "_peaks.annotatePeaks.txt",
            params.modules['plot_homer_annotatepeaks']
        )

        // Create channel: [ meta , [ peaks ] ]
        // Where meta = [ id:antibody, multiple_groups:true/false, replicates_exist:true/false ]
        MACS2_CALLPEAK
            .out
            .peak
            .map { meta, peak -> [ meta.antibody , meta.id.split('_')[0..-2].join('_'), peak ] }
            .groupTuple()
            .map {
                antibody, groups, peaks ->
                    [ antibody,
                      groups.groupBy().collectEntries { [(it.key) : it.value.size()] },
                      peaks ] }
            .map {
                antibody, groups, peaks ->
                    def meta = [:]
                    meta.id = antibody
                    meta.multiple_groups = groups.size() > 1
                    meta.replicates_exist = groups.max { groups.value }.value > 1
                    [ meta, peaks ] }
            .set { ch_antibody_peaks }

        params.modules['macs2_consensus'].publish_dir += "/$peakType/consensus"
        MACS2_CONSENSUS (
            ch_antibody_peaks,
            params.modules['macs2_consensus']
        )

        params.modules['homer_annotatepeaks_consensus'].publish_dir += "/$peakType/consensus"
        HOMER_ANNOTATEPEAKS_CONSENSUS (
            MACS2_CONSENSUS.out.bed,
            ch_fasta,
            ch_gtf,
            params.modules['homer_annotatepeaks_consensus']
        )

        //ch_software_versions = ch_software_versions.mix(SUBREAD_FEATURECOUNTS.out.version.first().ifEmpty(null))
    }


    /*
     * Pipeline reporting
     */
    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect(),
        params.modules['get_software_versions']
    )

    OUTPUT_DOCUMENTATION (
        ch_output_docs,
        ch_output_docs_images,
        params.modules['output_documentation']
    )

    /*
     * MultiQC
     */
    workflow_summary = Schema.params_mqc_summary(summary)
    ch_workflow_summary = Channel.value(workflow_summary)
    params.modules['multiqc'].publish_dir += "/$peakType"
    MULTIQC (
        ch_multiqc_config,
        ch_multiqc_custom_config.collect().ifEmpty([]),
        GET_SOFTWARE_VERSIONS.out.yaml.collect(),
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),

        FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
        FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),
        FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]),

        MAP_BWA_MEM.out.stats.collect{it[1]},
        MAP_BWA_MEM.out.flagstat.collect{it[1]},
        MAP_BWA_MEM.out.idxstats.collect{it[1]},

        MARK_DUPLICATES_PICARD.out.stats.collect{it[1]}.ifEmpty([]),
        MARK_DUPLICATES_PICARD.out.flagstat.collect{it[1]}.ifEmpty([]),
        MARK_DUPLICATES_PICARD.out.idxstats.collect{it[1]}.ifEmpty([]),
        MARK_DUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]),

        BAM_CLEAN.out.stats.collect{it[1]}.ifEmpty([]),
        BAM_CLEAN.out.flagstat.collect{it[1]}.ifEmpty([]),
        BAM_CLEAN.out.idxstats.collect{it[1]}.ifEmpty([]),
        PICARD_COLLECTMULTIPLEMETRICS.out.metrics.collect{it[1]}.ifEmpty([]),

        PRESEQ_LCEXTRAP.out.ccurve.collect{it[1]}.ifEmpty([]),
        DEEPTOOLS_PLOTPROFILE.out.table.collect{it[1]}.ifEmpty([]),
        DEEPTOOLS_PLOTFINGERPRINT.out.matrix.collect{it[1]}.ifEmpty([]),
        // path ('phantompeakqualtools/*') from ch_spp_out_mqc.collect().ifEmpty([])
        // path ('phantompeakqualtools/*') from ch_spp_csv_mqc.collect().ifEmpty([])

        // path ('macs/*') from ch_macs_mqc.collect().ifEmpty([])
        // path ('macs/*') from ch_macs_qc_mqc.collect().ifEmpty([])
        // path ('macs/consensus/*') from ch_macs_consensus_counts_mqc.collect().ifEmpty([])
        // path ('macs/consensus/*') from ch_macs_consensus_deseq_mqc.collect().ifEmpty([])

        params.modules['multiqc']
    )
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    def multiqc_report = []
    Completion.email(workflow, params, summary, run_name, baseDir, multiqc_report, log)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////

//
// /*
//  * STEP 7.2: Annotate consensus peaks with HOMER, and add annotation to boolean output file
//  */
// process CONSENSUS_PEAKS_ANNOTATE {
//     tag "${antibody}"
//     label 'process_medium'
//     publishDir "${params.outdir}/bwa/mergedLibrary/macs/${PEAK_TYPE}/consensus/${antibody}", mode: params.publish_dir_mode
//
//     when:
//     params.macs_gsize && (replicatesExist || multipleGroups) && !params.skip_consensus_peaks && !params.skip_peak_annotation
//
//     input:
//     tuple val(antibody), val(replicatesExist), val(multipleGroups), path(bed) from ch_macs_consensus_bed
//     path bool from ch_macs_consensus_bool
//     path fasta from ch_fasta
//     path gtf from ch_gtf
//
//     output:
//     path '*.annotatePeaks.txt'
//
//     script:
//     prefix = "${antibody}.consensus_peaks"
//     """
//     annotatePeaks.pl \\
//         $bed \\
//         $fasta \\
//         -gid \\
//         -gtf $gtf \\
//         -cpu $task.cpus \\
//         > ${prefix}.annotatePeaks.txt
//
//     cut -f2- ${prefix}.annotatePeaks.txt | awk 'NR==1; NR > 1 {print \$0 | "sort -T '.' -k1,1 -k2,2n"}' | cut -f6- > tmp.txt
//     paste $bool tmp.txt > ${prefix}.boolean.annotatePeaks.txt
//     """
// }
//
// // Get BAM and SAF files for each ip
// ch_group_bam_counts
//     .map { it -> [ it[3], [ it[0], it[1], it[2] ] ] }
//     .join(ch_rm_orphan_name_bam_counts)
//     .map { it -> [ it[1][0], it[1][1], it[1][2], it[2] ] }
//     .groupTuple()
//     .map { it -> [ it[0], it[1][0], it[2][0], it[3].flatten().sort() ] }
//     .join(ch_macs_consensus_saf)
//     .set { ch_group_bam_counts }
//
// /*
//  * STEP 7.3: Count reads in consensus peaks with featureCounts
//  */
// process CONSENSUS_PEAKS_COUNTS {
//     tag "${antibody}"
//     label 'process_medium'
//     publishDir "${params.outdir}/bwa/mergedLibrary/macs/${PEAK_TYPE}/consensus/${antibody}", mode: params.publish_dir_mode
//
//     when:
//     params.macs_gsize && (replicatesExist || multipleGroups) && !params.skip_consensus_peaks
//
//     input:
//     tuple val(antibody), val(replicatesExist), val(multipleGroups), path(bams), path(saf) from ch_group_bam_counts
//
//     output:
//     tuple val(antibody), val(replicatesExist), val(multipleGroups), path('*featureCounts.txt') into ch_macs_consensus_counts
//     path '*featureCounts.txt.summary' into ch_macs_consensus_counts_mqc
//
//     script:
//     prefix = "${antibody}.consensus_peaks"
//     bam_files = bams.findAll { it.toString().endsWith('.bam') }.sort()
//     pe_params = params.single_end ? '' : '-p --donotsort'
//     """
//     featureCounts \\
//         -F SAF \\
//         -O \\
//         --fracOverlap 0.2 \\
//         -T $task.cpus \\
//         $pe_params \\
//         -a $saf \\
//         -o ${prefix}.featureCounts.txt \\
//         ${bam_files.join(' ')}
//     """
// }
//
// /*
//  * STEP 7.4: Differential analysis with DESeq2
//  */
// process CONSENSUS_PEAKS_DESEQ2 {
//     tag "${antibody}"
//     label 'process_medium'
//     publishDir "${params.outdir}/bwa/mergedLibrary/macs/${PEAK_TYPE}/consensus/${antibody}/deseq2", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith('.igv.txt')) null
//                       else filename
//                 }
//
//     when:
//     params.macs_gsize && replicatesExist && multipleGroups && !params.skip_consensus_peaks && !params.skip_diff_analysis
//
//     input:
//     tuple val(antibody), val(replicatesExist), val(multipleGroups), path(counts) from ch_macs_consensus_counts
//     path deseq2_pca_header from ch_deseq2_pca_header
//     path deseq2_clustering_header from ch_deseq2_clustering_header
//
//     output:
//     path '*.tsv' into ch_macs_consensus_deseq_mqc
//     path '*igv.txt' into ch_macs_consensus_deseq_comp_igv
//     path '*.{RData,results.txt,pdf,log}'
//     path 'sizeFactors'
//     path '*vs*/*.{pdf,txt}'
//     path '*vs*/*.bed'
//
//     script:
//     prefix = "${antibody}.consensus_peaks"
//     bam_ext = params.single_end ? '.mLb.clN.sorted.bam' : '.mLb.clN.bam'
//     vst = params.deseq2_vst ? '--vst TRUE' : ''
//     """
//     featurecounts_deseq2.r \\
//         --featurecount_file $counts \\
//         --bam_suffix '$bam_ext' \\
//         --outdir ./ \\
//         --outprefix $prefix \\
//         --outsuffix '' \\
//         --cores $task.cpus \\
//         $vst
//
//     sed 's/deseq2_pca/deseq2_pca_${task.index}/g' <$deseq2_pca_header >tmp.txt
//     sed -i -e 's/DESeq2 /${antibody} DESeq2 /g' tmp.txt
//     cat tmp.txt ${prefix}.pca.vals.txt > ${prefix}.pca.vals_mqc.tsv
//
//     sed 's/deseq2_clustering/deseq2_clustering_${task.index}/g' <$deseq2_clustering_header >tmp.txt
//     sed -i -e 's/DESeq2 /${antibody} DESeq2 /g' tmp.txt
//     cat tmp.txt ${prefix}.sample.dists.txt > ${prefix}.sample.dists_mqc.tsv
//
//     find * -type f -name "*.FDR0.05.results.bed" -exec echo -e "bwa/mergedLibrary/macs/${PEAK_TYPE}/consensus/${antibody}/deseq2/"{}"\\t255,0,0" \\; > ${prefix}.igv.txt
//     """
// }
//
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                             IGV                                     -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
//
// /*
//  * STEP 8: Create IGV session file
//  */
// process IGV {
//     publishDir "${params.outdir}/igv/${PEAK_TYPE}", mode: params.publish_dir_mode
//
//     when:
//     !params.skip_igv
//
//     input:
//     path fasta from ch_fasta
//     path bigwigs from ch_bigwig_igv.collect().ifEmpty([])
//     path peaks from ch_macs_igv.collect().ifEmpty([])
//     path consensus_peaks from ch_macs_consensus_igv.collect().ifEmpty([])
//     path differential_peaks from ch_macs_consensus_deseq_comp_igv.collect().ifEmpty([])
//
//     output:
//     path '*.{txt,xml}'
//
//     script: // scripts are bundled with the pipeline in nf-core/chipseq/bin/
//     """
//     cat *.txt > igv_files.txt
//     igv_files_to_session.py igv_session.xml igv_files.txt ../../genome/${fasta.getName()} --path_prefix '../../'
//     """
// }

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        END OF PIPELINE                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
