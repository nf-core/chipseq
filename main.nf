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

nextflow.preview.dsl = 2

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

// Global variables
def PEAK_TYPE = params.narrow_peak ? 'narrowPeak' : 'broadPeak'

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
if (params.single_end) {
    ch_bamtools_filter_config = file(params.bamtools_filter_se_config, checkIfExists: true)
} else {
    ch_bamtools_filter_config = file(params.bamtools_filter_pe_config, checkIfExists: true)
}

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

include { GTF2BED               } from './modules/local/process/gtf2bed'
include { GET_CHROM_SIZES       } from './modules/local/process/get_chrom_sizes'
include { MAKE_GENOME_FILTER    } from './modules/local/process/make_genome_filter'
include { BEDTOOLS_GENOMECOV    } from './modules/local/process/bedtools_genomecov'
include { OUTPUT_DOCUMENTATION  } from './modules/local/process/output_documentation'
include { GET_SOFTWARE_VERSIONS } from './modules/local/process/get_software_versions'

include { CHECK_INPUT           } from './modules/local/subworkflow/check_input'
include { CLEAN_BAM             } from './modules/local/subworkflow/clean_bam'

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
//include { MACSC2_CALLPEAK             } from './modules/nf-core/software/macs2_callpeak'
//include { HOMER_ANNOTATEPEAKS         } from './modules/nf-core/software/homer_annotatepeaks'
//include { SUBREAD_FEATURECOUNTS       } from './modules/nf-core/software/subread_featurecounts'
include { MULTIQC                       } from './modules/nf-core/software/multiqc'

include { QC_TRIM                       } from './modules/nf-core/subworkflow/qc_trim'
include { BAM_STATS                     } from './modules/nf-core/subworkflow/bam_stats'
include { SORT_BAM                      } from './modules/nf-core/subworkflow/sort_bam'
include { MAP_READS                     } from './modules/nf-core/subworkflow/map_reads'
include { MARK_DUPLICATES               } from './modules/nf-core/subworkflow/mark_duplicates'

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {

    /*
     * Read in samplesheet, validate and stage input files
     */
    CHECK_INPUT (
        ch_input,
        params.seq_center,
        params.modules['check_samplesheet']
    )

    /*
     * Prepare genome files
     */
    ch_index = params.bwa_index ? Channel.value(file(params.bwa_index)) : BWA_INDEX ( ch_fasta, params.modules['bwa_index'] ).index

    if (makeBED) { ch_gene_bed = GTF2BED ( ch_gtf, params.modules['gtf2bed'] ) }

    MAKE_GENOME_FILTER (
        GET_CHROM_SIZES (
            ch_fasta,
            params.modules['get_chrom_sizes']
        ).sizes,
        ch_blacklist.ifEmpty([]),
        params.modules['make_genome_filter']
    )

    /*
     * Read QC & trimming
     */
    nextseq = params.trim_nextseq > 0 ? " --nextseq ${params.trim_nextseq}" : ''
    params.modules['trimgalore'].args += nextseq
    QC_TRIM (
        CHECK_INPUT.out.reads,
        params.skip_fastqc,
        params.skip_trimming,
        params.modules['fastqc'],
        params.modules['trimgalore']
    )

    /*
     * Map reads & BAM QC
     */
    score = params.bwa_min_score ? " -T ${params.bwa_min_score}" : ''
    params.modules['bwa_mem'].args += score
    MAP_READS (
        QC_TRIM.out.reads,
        ch_index,
        ch_fasta,
        params.modules['bwa_mem'],
        params.modules['samtools_sort_lib']
    )

    /*
     * Merge resequenced BAM files
     */
    MAP_READS
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

    /*
     * Mark duplicates & filter BAM files
     */
    MARK_DUPLICATES (
        PICARD_MERGESAMFILES.out.bam,
        params.modules['picard_markduplicates'],
        params.modules['samtools_sort_merged_lib']
    )

    // Fix getting name sorted BAM here for PE/SE
    CLEAN_BAM (
        MARK_DUPLICATES.out.bam.join(MARK_DUPLICATES.out.bai, by: [0]),
        MAKE_GENOME_FILTER.out.collect(),
        ch_bamtools_filter_config,
        params.modules['filter_bam'],
        params.modules['remove_bam_orphans'],
        params.modules['samtools_sort_filter']
    )

    /*
     * Post alignment QC
     */
    PICARD_COLLECTMULTIPLEMETRICS (
        CLEAN_BAM.out.bam,
        ch_fasta,
        params.modules['picard_collectmultiplemetrics']
    )

    PRESEQ_LCEXTRAP (
        CLEAN_BAM.out.bam,
        params.modules['preseq_lcextrap']
    )

    PHANTOMPEAKQUALTOOLS (
        CLEAN_BAM.out.bam,
        params.modules['phantompeakqualtools']
    )

    /*
     * Coverage tracks
     */
    BEDTOOLS_GENOMECOV (
        CLEAN_BAM.out.bam.join(CLEAN_BAM.out.flagstat, by: [0]),
        params.modules['bedtools_genomecov']
    )

    UCSC_BEDRAPHTOBIGWIG (
        BEDTOOLS_GENOMECOV.out.bedgraph,
        GET_CHROM_SIZES.out.sizes,
        params.modules['ucsc_bedgraphtobigwig']
    )

    /*
     * Coverage plots
     */
    DEEPTOOLS_COMPUTEMATRIX (
        UCSC_BEDRAPHTOBIGWIG.out.bigwig,
        ch_gene_bed,
        params.modules['deeptools_computematrix']
    )

    DEEPTOOLS_PLOTPROFILE (
        DEEPTOOLS_COMPUTEMATRIX.out.matrix,
        params.modules['deeptools_plotprofile']
    )

    DEEPTOOLS_PLOTHEATMAP (
        DEEPTOOLS_COMPUTEMATRIX.out.matrix,
        params.modules['deeptools_plotheatmap']
    )

    // Join control BAM here too to generate plots with IP and CONTROL together
    params.modules['deeptools_plotfingerprint'].args += " --numberOfSamples $params.fingerprint_bins"
    DEEPTOOLS_PLOTFINGERPRINT (
        CLEAN_BAM.out.bam.join(CLEAN_BAM.out.bai, by: [0]),
        params.modules['deeptools_plotfingerprint']
    )

    /*
     * Pipeline reporting
     */
    GET_SOFTWARE_VERSIONS (
        params.modules['get_software_versions']
    )

    OUTPUT_DOCUMENTATION (
        ch_output_docs,
        ch_output_docs_images,
        params.modules['output_documentation']
    )

    // MULTIQC(
    //     summary,
    //     FASTQC.out,
    //     ch_multiqc_config
    // )
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

// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                 MERGE LIBRARY PEAK ANALYSIS                         -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
//
// // Create channel linking IP bams with control bams
// ch_rm_orphan_bam_macs_1
//     .combine(ch_rm_orphan_bam_macs_2)
//     .set { ch_rm_orphan_bam_macs_1 }
//
// ch_design_controls_csv
//     .combine(ch_rm_orphan_bam_macs_1)
//     .filter { it[0] == it[5] && it[1] == it[7] }
//     .join(ch_rm_orphan_flagstat_macs)
//     .map { it ->  it[2..-1] }
//     .into { ch_group_bam_macs;
//             ch_group_bam_plotfingerprint;
//             ch_group_bam_counts }
//
// /*
//  * STEP 6.2: Call peaks with MACS2 and calculate FRiP score
//  */
// process MACS2 {
//     tag "${ip} vs ${control}"
//     label 'process_medium'
//     publishDir "${params.outdir}/bwa/mergedLibrary/macs/${PEAK_TYPE}", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith('.tsv')) "qc/$filename"
//                       else if (filename.endsWith('.igv.txt')) null
//                       else filename
//                 }
//
//     when:
//     params.macs_gsize
//
//     input:
//     tuple val(antibody), val(replicatesExist), val(multipleGroups), val(ip), path(ipbam), val(control), path(controlbam), path(ipflagstat) from ch_group_bam_macs
//     path peak_count_header from ch_peak_count_header
//     path frip_score_header from ch_frip_score_header
//
//     output:
//     tuple val(antibody), val(replicatesExist), val(multipleGroups), val(ip), val(control), path("*.$PEAK_TYPE") into ch_macs_homer,
//                                                                                                                      ch_macs_qc,
//                                                                                                                      ch_macs_consensus
//     path '*igv.txt' into ch_macs_igv
//     path '*_mqc.tsv' into ch_macs_mqc
//     path '*.{bed,xls,gappedPeak,bdg}'
//
//     script:
//     broad = params.narrow_peak ? '' : "--broad --broad-cutoff ${params.broad_cutoff}"
//     format = params.single_end ? 'BAM' : 'BAMPE'
//     pileup = params.save_macs_pileup ? '-B --SPMR' : ''
//     fdr = params.macs_fdr ? "--qvalue ${params.macs_fdr}" : ''
//     pvalue = params.macs_pvalue ? "--pvalue ${params.macs_pvalue}" : ''
//     """
//     macs2 callpeak \\
//         -t ${ipbam[0]} \\
//         -c ${controlbam[0]} \\
//         $broad \\
//         -f $format \\
//         -g $params.macs_gsize \\
//         -n $ip \\
//         $pileup \\
//         $fdr \\
//         $pvalue \\
//         --keep-dup all
//
//     cat ${ip}_peaks.${PEAK_TYPE} | wc -l | awk -v OFS='\t' '{ print "${ip}", \$1 }' | cat $peak_count_header - > ${ip}_peaks.count_mqc.tsv
//
//     READS_IN_PEAKS=\$(intersectBed -a ${ipbam[0]} -b ${ip}_peaks.${PEAK_TYPE} -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
//     grep 'mapped (' $ipflagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${ip}", a/\$1}' | cat $frip_score_header - > ${ip}_peaks.FRiP_mqc.tsv
//
//     find * -type f -name "*.${PEAK_TYPE}" -exec echo -e "bwa/mergedLibrary/macs/${PEAK_TYPE}/"{}"\\t0,0,178" \\; > ${ip}_peaks.igv.txt
//     """
// }
//
// /*
//  * STEP 6.3: Annotate peaks with HOMER
//  */
// process MACS2_ANNOTATE {
//     tag "${ip} vs ${control}"
//     label 'process_medium'
//     publishDir "${params.outdir}/bwa/mergedLibrary/macs/${PEAK_TYPE}", mode: params.publish_dir_mode
//
//     when:
//     params.macs_gsize && !params.skip_peak_annotation
//
//     input:
//     tuple val(antibody), val(replicatesExist), val(multipleGroups), val(ip), val(control), path(peak) from ch_macs_homer
//     path fasta from ch_fasta
//     path gtf from ch_gtf
//
//     output:
//     path '*.txt' into ch_macs_annotate
//
//     script:
//     """
//     annotatePeaks.pl \\
//         $peak \\
//         $fasta \\
//         -gid \\
//         -gtf $gtf \\
//         -cpu $task.cpus \\
//         > ${ip}_peaks.annotatePeaks.txt
//     """
// }
//
// /*
//  * STEP 6.4: Aggregated QC plots for peaks, FRiP and peak-to-gene annotation
//  */
// process MACS2_QC {
//     label 'process_medium'
//     publishDir "${params.outdir}/bwa/mergedLibrary/macs/${PEAK_TYPE}/qc", mode: params.publish_dir_mode
//
//     when:
//     params.macs_gsize && !params.skip_peak_annotation && !params.skip_peak_qc
//
//     input:
//     path peaks from ch_macs_qc.collect{ it[-1] }
//     path annos from ch_macs_annotate.collect()
//     path peak_annotation_header from ch_peak_annotation_header
//
//     output:
//     path '*.tsv' into ch_macs_qc_mqc
//     path '*.{txt,pdf}'
//
//     script:  // This script is bundled with the pipeline, in nf-core/chipseq/bin/
//     """
//     plot_macs_qc.r \\
//         -i ${peaks.join(',')} \\
//         -s ${peaks.join(',').replaceAll("_peaks.${PEAK_TYPE}","")} \\
//         -o ./ \\
//         -p macs_peak
//
//     plot_homer_annotatepeaks.r \\
//         -i ${annos.join(',')} \\
//         -s ${annos.join(',').replaceAll("_peaks.annotatePeaks.txt","")} \\
//         -o ./ \\
//         -p macs_annotatePeaks
//
//     cat $peak_annotation_header macs_annotatePeaks.summary.txt > macs_annotatePeaks.summary_mqc.tsv
//     """
// }
//
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                 CONSENSUS PEAKS ANALYSIS                            -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
//
// // Group by ip from this point and carry forward boolean variables
// ch_macs_consensus
//     .map { it ->  [ it[0], it[1], it[2], it[-1] ] }
//     .groupTuple()
//     .map { it ->  [ it[0], it[1][0], it[2][0], it[3].sort() ] }
//     .set { ch_macs_consensus }
//
// /*
//  * STEP 7.1: Consensus peaks across samples, create boolean filtering file, SAF file for featureCounts and UpSetR plot for intersection
//  */
// process CONSENSUS_PEAKS {
//     tag "${antibody}"
//     label 'process_long'
//     publishDir "${params.outdir}/bwa/mergedLibrary/macs/${PEAK_TYPE}/consensus/${antibody}", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith('.igv.txt')) null
//                       else filename
//                 }
//
//     when:
//     params.macs_gsize && (replicatesExist || multipleGroups) && !params.skip_consensus_peaks
//
//     input:
//     tuple val(antibody), val(replicatesExist), val(multipleGroups), path(peaks) from ch_macs_consensus
//
//     output:
//     tuple val(antibody), val(replicatesExist), val(multipleGroups), path('*.bed') into ch_macs_consensus_bed
//     tuple val(antibody), path('*.saf') into ch_macs_consensus_saf
//     path '*.boolean.txt' into ch_macs_consensus_bool
//     path '*igv.txt' into ch_macs_consensus_igv
//     path '*.intersect.{txt,plot.pdf}'
//
//     script: // scripts are bundled with the pipeline in nf-core/chipseq/bin/
//     prefix = "${antibody}.consensus_peaks"
//     mergecols = params.narrow_peak ? (2..10).join(',') : (2..9).join(',')
//     collapsecols = params.narrow_peak ? (['collapse']*9).join(',') : (['collapse']*8).join(',')
//     expandparam = params.narrow_peak ? '--is_narrow_peak' : ''
//     """
//     sort -T '.' -k1,1 -k2,2n ${peaks.collect{it.toString()}.sort().join(' ')} \\
//         | mergeBed -c $mergecols -o $collapsecols > ${prefix}.txt
//
//     macs2_merged_expand.py ${prefix}.txt \\
//         ${peaks.collect{it.toString()}.sort().join(',').replaceAll("_peaks.${PEAK_TYPE}","")} \\
//         ${prefix}.boolean.txt \\
//         --min_replicates $params.min_reps_consensus \\
//         $expandparam
//
//     awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2, \$3, \$4, "0", "+" }' ${prefix}.boolean.txt > ${prefix}.bed
//
//     echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${prefix}.saf
//     awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$4, \$1, \$2, \$3,  "+" }' ${prefix}.boolean.txt >> ${prefix}.saf
//
//     plot_peak_intersect.r -i ${prefix}.boolean.intersect.txt -o ${prefix}.boolean.intersect.plot.pdf
//
//     find * -type f -name "${prefix}.bed" -exec echo -e "bwa/mergedLibrary/macs/${PEAK_TYPE}/consensus/${antibody}/"{}"\\t0,0,0" \\; > ${prefix}.bed.igv.txt
//     """
// }
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
//
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                          MULTIQC                                    -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
//
// Channel.from(summary.collect{ [it.key, it.value] })
//     .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
//     .reduce { a, b -> return [a, b].join("\n            ") }
//     .map { x -> """
//     id: 'nf-core-chipseq-summary'
//     description: " - this information is collected when the pipeline is started."
//     section_name: 'nf-core/chipseq Workflow Summary'
//     section_href: 'https://github.com/nf-core/chipseq'
//     plot_type: 'html'
//     data: |
//         <dl class=\"dl-horizontal\">
//             $x
//         </dl>
//     """.stripIndent() }
//     .set { ch_workflow_summary }
//
// /*
//  * STEP 9: MultiQC
//  */
// process MULTIQC {
//     publishDir "${params.outdir}/multiqc/${PEAK_TYPE}", mode: params.publish_dir_mode
//
//     when:
//     !params.skip_multiqc
//
//     input:
//     path (multiqc_config) from ch_multiqc_config
//     path (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
//
//     path ('software_versions/*') from ch_software_versions_mqc.collect()
//     path workflow_summary from ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
//
//     path ('fastqc/*') from ch_fastqc_reports_mqc.collect().ifEmpty([])
//     path ('trimgalore/*') from ch_trimgalore_results_mqc.collect().ifEmpty([])
//     path ('trimgalore/fastqc/*') from ch_trimgalore_fastqc_reports_mqc.collect().ifEmpty([])
//
//     path ('alignment/library/*') from ch_sort_bam_flagstat_mqc.collect()
//     path ('alignment/mergedLibrary/*') from ch_merge_bam_stats_mqc.collect()
//     path ('alignment/mergedLibrary/*') from ch_rm_orphan_flagstat_mqc.collect{it[1]}
//     path ('alignment/mergedLibrary/*') from ch_rm_orphan_stats_mqc.collect()
//     path ('alignment/mergedLibrary/picard_metrics/*') from ch_merge_bam_metrics_mqc.collect()
//     path ('alignment/mergedLibrary/picard_metrics/*') from ch_collectmetrics_mqc.collect()
//
//     path ('macs/*') from ch_macs_mqc.collect().ifEmpty([])
//     path ('macs/*') from ch_macs_qc_mqc.collect().ifEmpty([])
//     path ('macs/consensus/*') from ch_macs_consensus_counts_mqc.collect().ifEmpty([])
//     path ('macs/consensus/*') from ch_macs_consensus_deseq_mqc.collect().ifEmpty([])
//
//     path ('preseq/*') from ch_preseq_mqc.collect().ifEmpty([])
//     path ('deeptools/*') from ch_plotfingerprint_mqc.collect().ifEmpty([])
//     path ('deeptools/*') from ch_plotprofile_mqc.collect().ifEmpty([])
//     path ('phantompeakqualtools/*') from ch_spp_out_mqc.collect().ifEmpty([])
//     path ('phantompeakqualtools/*') from ch_spp_csv_mqc.collect().ifEmpty([])
//
//     output:
//     path '*multiqc_report.html' into ch_multiqc_report
//     path '*_data'
//
//     script:
//     rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
//     rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
//     custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
//     """
//     multiqc . -f $rtitle $rfilename $custom_config_file
//     """
// }
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        END OF PIPELINE                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
