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
    lastPath = params.fasta.lastIndexOf(File.separator)
    params.bwa_base = params.fasta.substring(lastPath+1)
    ch_fasta = file(params.fasta, checkIfExists: true)
} else {
    exit 1, 'Fasta file not specified!'
}

if (params.bwa_index) {
    lastPath = params.bwa_index.lastIndexOf(File.separator)
    bwa_dir  = params.bwa_index.substring(0,lastPath+1)
    params.bwa_base = params.bwa_index.substring(lastPath+1)
    Channel
        .fromPath(bwa_dir, checkIfExists: true)
        .set { ch_index }
}

// Save AWS IGenomes file containing annotation version
if (params.anno_readme && file(params.anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(params.anno_readme).copyTo("${params.outdir}/genome/")
}

// If --gtf is supplied along with --genome
// Make gene bed from supplied --gtf instead of using iGenomes one automatically
def MakeBED = false
if (!params.gene_bed) {
    MakeBED = true
} else if (params.genome && params.gtf) {
    if (params.genomes[ params.genome ].gtf != params.gtf) {
        MakeBED = true
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

/*
 * Include local pipeline modules
 */
include { CHECK_SAMPLESHEET;
          get_samplesheet_paths;
          get_samplesheet_design } from './modules/local/check_samplesheet'
include { GTF2BED } from './modules/local/gtf2bed'
include { GET_CHROM_SIZES } from './modules/local/get_chrom_sizes'
include { MAKE_GENOME_FILTER } from './modules/local/make_genome_filter'
include { FILTER_BAM } from './modules/local/filter_bam'
include { OUTPUT_DOCUMENTATION } from './modules/local/output_documentation'
include { GET_SOFTWARE_VERSIONS } from './modules/local/get_software_versions'

/*
 * Include nf-core modules
 */
include { FASTQC } from './modules/nf-core/fastqc' addParams(fastqc_args : "--quiet")
include { BWA_INDEX } from './modules/nf-core/bwa_index'
include { TRIMGALORE } from './modules/nf-core/trimgalore'
include { BWA_MEM } from './modules/nf-core/bwa_mem'
include { SAMTOOLS_SORT } from './modules/nf-core/samtools_sort'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools_index'
include { SAMTOOLS_STATS } from './modules/nf-core/samtools_stats'
include { PICARD_MERGESAMFILES } from './modules/nf-core/picard_mergesamfiles'
include { PICARD_MARKDUPLICATES } from './modules/nf-core/picard_markduplicates'
//include { PICARD_COLLECTMULTIPLEMETRICS } from './modules/nf-core/picard_collectmultiplemetrics'
include { MULTIQC } from './modules/nf-core/multiqc'

/*
 * Check input samplesheet and get read channels
 */
workflow CHECK_INPUT {
    take:
    ch_input // file: /path/to/samplesheet.csv

    main:
    CHECK_SAMPLESHEET(ch_input)
        .reads
        .splitCsv(header:true, sep:',')
        .map { get_samplesheet_paths(it, params.single_end) }
        .set { ch_reads }

    CHECK_SAMPLESHEET
        .out
        .controls
        .splitCsv(header:true, sep:',')
        .map { get_samplesheet_design(it, params.single_end) }
        .set { ch_design }

    emit:
    reads = ch_reads
    design = ch_design
}

/*
 * Read QC and trimming
 */
workflow QC_TRIM {
    take:
    ch_reads        // channel: [ val(metadata), [ reads ] ]
    skip_fastqc     // boolean: true/false
    skip_trimming   // boolean: true/false
    fastqc_opts     //     map: options for FastQC module
    trimgalore_opts //     map: options for TrimGalore! module

    main:
    fastqc_html = Channel.empty()
    fastqc_zip = Channel.empty()
    fastqc_version = Channel.empty()
    if (!skip_fastqc) {
        FASTQC(ch_reads, fastqc_opts).html.set { fastqc_html }
        fastqc_zip = FASTQC.out.zip
        fastqc_version = FASTQC.out.version
    }

    // ch_trim_reads = ch_reads
    // trim_log = Channel.empty()
    // trim_fastqc = Channel.empty()
    // if (!skip_trimming) {
    //     TRIMGALORE(ch_reads, trimgalore_opts).reads.set { ch_trim_reads }
    //     trim_log = TRIMGALORE.out.log
    //     trim_fastqc = TRIMGALORE.out.fastqc
    // }

    emit:
    fastqc_html = fastqc_html
    fastqc_zip = fastqc_zip
    fastqc_version = fastqc_version
    // reads = ch_trim_reads
    // trim_log = trim_log
    // trim_fastqc = trim_fastqc
}

/*
 * Sort, index BAM file and run samtools stats, flagstat and idxstats
 */
workflow SORT_BAM {
    take:
    ch_bam // channel: [ val(name), val(single_end), bam ]

    main:
    SAMTOOLS_SORT(ch_bam) | SAMTOOLS_INDEX
    SAMTOOLS_STATS(SAMTOOLS_SORT.out.join(SAMTOOLS_INDEX.out, by: [0,1]))

    emit:
    bam = SAMTOOLS_SORT.out
    bai = SAMTOOLS_INDEX.out
    stats = SAMTOOLS_STATS.out.stats
    flagstat = SAMTOOLS_STATS.out.flagstat
    idxstats = SAMTOOLS_STATS.out.idxstats
}

/*
 * Map reads, sort, index BAM file and run samtools stats, flagstat and idxstats
 */
workflow ALIGN {
    take:
    ch_reads // channel: [ val(name), val(single_end), [ reads ] ]
    ch_index // channel: [ /path/to/index ]

    main:
    BWA_MEM(ch_reads, ch_index.collect())
    SORT_BAM(BWA_MEM.out)

    emit:
    bam = SORT_BAM.out.bam
    bai = SORT_BAM.out.bai
    stats = SORT_BAM.out.stats
    flagstat = SORT_BAM.out.flagstat
    idxstats = SORT_BAM.out.idxstats
}

/*
 * Picard MarkDuplicates, sort, index BAM file and run samtools stats, flagstat and idxstats
 */
workflow MARKDUP {
    take:
    ch_bam // channel: [ val(name), val(single_end), bam ]

    main:
    PICARD_MARKDUPLICATES(ch_bam)
    SAMTOOLS_INDEX(PICARD_MARKDUPLICATES.out.bam)
    SAMTOOLS_STATS(PICARD_MARKDUPLICATES.out.bam.join(SAMTOOLS_INDEX.out, by: [0,1]))

    emit:
    bam = PICARD_MARKDUPLICATES.out.bam
    metrics = PICARD_MARKDUPLICATES.out.metrics
    bai = SAMTOOLS_INDEX.out
    stats = SAMTOOLS_STATS.out.stats
    flagstat = SAMTOOLS_STATS.out.flagstat
    idxstats = SAMTOOLS_STATS.out.idxstats
}

/*
 * Filter BAM file
 */
workflow FILTER_BAM {
    take:
    ch_bam // channel: [ val(name), val(single_end), bam ]
    ch_bai // channel: [ val(name), val(single_end), bai ]
    ch_bed // channel: [ bed ]
    config //    file: BAMtools filter JSON config file

    main:
    FILTER_BAM(ch_bam,
               ch_bai,
               ch_bed,
               config)                  // Fix getting name sorted BAM here for PE/SE
    SAMTOOLS_INDEX(FILTER_BAM.out.bam)
    SAMTOOLS_STATS(FILTER_BAM.out.bam.join(SAMTOOLS_INDEX.out, by: [0,1]))

    emit:
    bam = FILTER_BAM.out.bam
    bai = SAMTOOLS_INDEX.out
    stats = SAMTOOLS_STATS.out.stats
    flagstat = SAMTOOLS_STATS.out.flagstat
    idxstats = SAMTOOLS_STATS.out.idxstats
}

workflow {

    // READ IN SAMPLESHEET, VALIDATE AND STAGE INPUT FILES
    CHECK_INPUT(ch_input)

    // // PREPARE GENOME FILES
    // if (!params.bwa_index) { BWA_INDEX(ch_fasta).set { ch_index } }
    // if (MakeBED) { GTF2BED(ch_gtf).set { ch_gene_bed } }
    // MAKE_GENOME_FILTER(GET_CHROM_SIZES(ch_fasta).sizes, ch_blacklist.ifEmpty([]))
    //
    // READ QC & TRIMMING
    def fastqc_opts = [:]
    fastqc_opts.args = "--quiet"
    fastqc_opts.suffix = "test"
    fastqc_opts.publish_dir = "fastqc_test"
    fastqc_opts.publish_results = "all"
    def trimgalore_opts = [:]
    QC_TRIM(CHECK_INPUT.out.reads, params.skip_fastqc, params.skip_trimming, fastqc_opts, trimgalore_opts)
    //
    // // MAP READS & BAM QC
    // ALIGN(QC_TRIM.out.reads, ch_index.collect())
    //
    // // MERGE RESEQUENCED BAM FILES
    // ALIGN
    //     .out
    //     .bam
    //     .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1], it[2] ] }
    //     .groupTuple(by: [0, 1])
    //     .map { it ->  [ it[0], it[1], it[2].flatten() ] }
    //     .set { ch_sort_bam }
    // PICARD_MERGESAMFILES(ch_sort_bam) | MARKDUP
    //
    // // FILTER BAM FILES
    // // FILTER_BAM(MARK_DUPS.out.bam,
    // //            MAKE_GENOME_FILTER.out.collect(),
    // //            ch_bamtools_filter_config)           // Fix getting name sorted BAM here for PE/SE
    //
    // // PIPELINE TEMPLATE REPORTING
    // GET_SOFTWARE_VERSIONS()
    // OUTPUT_DOCUMENTATION(ch_output_docs,ch_output_docs_images)
    //
    // // MULTIQC(
    // //     summary,
    // //     FASTQC.out,
    // //     ch_multiqc_config
    // // )
}

/*
 * Send completion email
 */
workflow.onComplete {
    def multiqc_report = []
    Completion.email(workflow, params, summary, run_name, baseDir, multiqc_report, log)
    Completion.summary(workflow, params, log)
}

// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                    MERGE LIBRARY BAM                                -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
//
// /*
//  * STEP 4.3: Remove orphan reads from paired-end BAM file
//  */
// if (params.single_end) {
//     ch_filter_bam
//         .into { ch_rm_orphan_bam_metrics;
//                 ch_rm_orphan_bam_bigwig;
//                 ch_rm_orphan_bam_macs_1;
//                 ch_rm_orphan_bam_macs_2;
//                 ch_rm_orphan_bam_phantompeakqualtools;
//                 ch_rm_orphan_name_bam_counts }
//
//     ch_filter_bam_flagstat
//         .into { ch_rm_orphan_flagstat_bigwig;
//                 ch_rm_orphan_flagstat_macs;
//                 ch_rm_orphan_flagstat_mqc }
//
//     ch_filter_bam_stats_mqc
//         .set { ch_rm_orphan_stats_mqc }
// } else {
//     process MERGED_BAM_REMOVE_ORPHAN {
//         tag "$name"
//         label 'process_medium'
//         publishDir path: "${params.outdir}/bwa/mergedLibrary", mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                           if (filename.endsWith('.flagstat')) "samtools_stats/$filename"
//                           else if (filename.endsWith('.idxstats')) "samtools_stats/$filename"
//                           else if (filename.endsWith('.stats')) "samtools_stats/$filename"
//                           else if (filename.endsWith('.sorted.bam')) filename
//                           else if (filename.endsWith('.sorted.bam.bai')) filename
//                           else null
//                     }
//
//         input:
//         tuple val(name), path(bam) from ch_filter_bam
//
//         output:
//         tuple val(name), path('*.sorted.{bam,bam.bai}') into ch_rm_orphan_bam_metrics,
//                                                              ch_rm_orphan_bam_bigwig,
//                                                              ch_rm_orphan_bam_macs_1,
//                                                              ch_rm_orphan_bam_macs_2,
//                                                              ch_rm_orphan_bam_phantompeakqualtools
//         tuple val(name), path("${prefix}.bam") into ch_rm_orphan_name_bam_counts
//         tuple val(name), path('*.flagstat') into ch_rm_orphan_flagstat_bigwig,
//                                                  ch_rm_orphan_flagstat_macs,
//                                                  ch_rm_orphan_flagstat_mqc
//         path '*.{idxstats,stats}' into ch_rm_orphan_stats_mqc
//
//         script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
//         prefix = "${name}.mLb.clN"
//         """
//         bampe_rm_orphan.py ${bam[0]} ${prefix}.bam --only_fr_pairs
//
//         samtools sort -@ $task.cpus -o ${prefix}.sorted.bam -T $prefix ${prefix}.bam
//         samtools index ${prefix}.sorted.bam
//         samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
//         samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
//         samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
//         """
//     }
// }
//
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                 MERGE LIBRARY BAM POST-ANALYSIS                     -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
//
// /*
//  * STEP 5.1: Preseq analysis after merging libraries and before filtering
//  */
// process PRESEQ {
//     tag "$name"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/bwa/mergedLibrary/preseq", mode: params.publish_dir_mode
//
//     when:
//     !params.skip_preseq
//
//     input:
//     tuple val(name), path(bam) from ch_merge_bam_preseq
//
//     output:
//     path '*.ccurve.txt' into ch_preseq_mqc
//     path '*.log'
//
//     script:
//     pe = params.single_end ? '' : '-pe'
//     """
//     preseq lc_extrap \\
//         -output ${name}.ccurve.txt \\
//         -verbose \\
//         -bam \\
//         $pe \\
//         -seed 1 \\
//         ${bam[0]}
//     cp .command.err ${name}.command.log
//     """
// }
//
// /*
//  * STEP 5.2: Picard CollectMultipleMetrics after merging libraries and filtering
//  */
// process PICARD_METRICS {
//     tag "$name"
//     label 'process_medium'
//     publishDir path: "${params.outdir}/bwa/mergedLibrary", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith('_metrics')) "picard_metrics/$filename"
//                       else if (filename.endsWith('.pdf')) "picard_metrics/pdf/$filename"
//                       else null
//                 }
//
//     when:
//     !params.skip_picard_metrics
//
//     input:
//     tuple val(name), path(bam) from ch_rm_orphan_bam_metrics
//     path fasta from ch_fasta
//
//     output:
//     path '*_metrics' into ch_collectmetrics_mqc
//     path '*.pdf'
//
//     script:
//     prefix = "${name}.mLb.clN"
//     def avail_mem = 3
//     if (!task.memory) {
//         log.info '[Picard CollectMultipleMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
//     } else {
//         avail_mem = task.memory.toGiga()
//     }
//     """
//     picard -Xmx${avail_mem}g CollectMultipleMetrics \\
//         INPUT=${bam[0]} \\
//         OUTPUT=${prefix}.CollectMultipleMetrics \\
//         REFERENCE_SEQUENCE=$fasta \\
//         VALIDATION_STRINGENCY=LENIENT \\
//         TMP_DIR=tmp
//     """
// }
//
// /*
//  * STEP 5.3: Read depth normalised bigWig
//  */
// process BIGWIG {
//     tag "$name"
//     label 'process_medium'
//     publishDir "${params.outdir}/bwa/mergedLibrary/bigwig", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith('scale_factor.txt')) "scale/$filename"
//                       else if (filename.endsWith('.bigWig')) filename
//                       else null
//                 }
//
//     input:
//     tuple val(name), path(bam), path(flagstat) from ch_rm_orphan_bam_bigwig.join(ch_rm_orphan_flagstat_bigwig, by: [0])
//     path sizes from ch_genome_sizes_bigwig.collect()
//
//     output:
//     tuple val(name), path('*.bigWig') into ch_bigwig_plotprofile
//     path '*igv.txt' into ch_bigwig_igv
//     path '*scale_factor.txt'
//
//     script:
//     pe_fragment = params.single_end ? '' : '-pc'
//     extend = (params.single_end && params.fragment_size > 0) ? "-fs ${params.fragment_size}" : ''
//     """
//     SCALE_FACTOR=\$(grep 'mapped (' $flagstat | awk '{print 1000000/\$1}')
//     echo \$SCALE_FACTOR > ${name}.scale_factor.txt
//     genomeCoverageBed -ibam ${bam[0]} -bg -scale \$SCALE_FACTOR $pe_fragment $extend | sort -T '.' -k1,1 -k2,2n >  ${name}.bedGraph
//
//     bedGraphToBigWig ${name}.bedGraph $sizes ${name}.bigWig
//
//     find * -type f -name "*.bigWig" -exec echo -e "bwa/mergedLibrary/bigwig/"{}"\\t0,0,178" \\; > ${name}.bigWig.igv.txt
//     """
// }
//
// /*
//  * STEP 5.4: Generate gene body coverage plot with deepTools plotProfile and plotHeatmap
//  */
// process PLOTPROFILE {
//     tag "$name"
//     label 'process_high'
//     publishDir "${params.outdir}/bwa/mergedLibrary/deepTools/plotProfile", mode: params.publish_dir_mode
//
//     when:
//     !params.skip_plot_profile
//
//     input:
//     tuple val(name), path(bigwig) from ch_bigwig_plotprofile
//     path bed from ch_gene_bed
//
//     output:
//     path '*.plotProfile.tab' into ch_plotprofile_mqc
//     path '*.{gz,pdf,mat.tab}'
//
//     script:
//     """
//     computeMatrix scale-regions \\
//         --regionsFileName $bed \\
//         --scoreFileName $bigwig \\
//         --outFileName ${name}.computeMatrix.mat.gz \\
//         --outFileNameMatrix ${name}.computeMatrix.vals.mat.tab \\
//         --regionBodyLength 1000 \\
//         --beforeRegionStartLength 3000 \\
//         --afterRegionStartLength 3000 \\
//         --skipZeros \\
//         --smartLabels \\
//         --numberOfProcessors $task.cpus
//
//     plotProfile --matrixFile ${name}.computeMatrix.mat.gz \\
//         --outFileName ${name}.plotProfile.pdf \\
//         --outFileNameData ${name}.plotProfile.tab
//
//     plotHeatmap --matrixFile ${name}.computeMatrix.mat.gz \\
//         --outFileName ${name}.plotHeatmap.pdf \\
//         --outFileNameMatrix ${name}.plotHeatmap.mat.tab
//     """
// }
//
// /*
//  * STEP 5.5: Phantompeakqualtools
//  */
// process PHANTOMPEAKQUALTOOLS {
//     tag "$name"
//     label 'process_medium'
//     publishDir "${params.outdir}/bwa/mergedLibrary/phantompeakqualtools", mode: params.publish_dir_mode
//
//     when:
//     !params.skip_spp
//
//     input:
//     tuple val(name), path(bam) from ch_rm_orphan_bam_phantompeakqualtools
//     path spp_correlation_header from ch_spp_correlation_header
//     path spp_nsc_header from ch_spp_nsc_header
//     path spp_rsc_header from ch_spp_rsc_header
//
//     output:
//     path '*.spp.out' into ch_spp_out_mqc
//     path '*_mqc.tsv' into ch_spp_csv_mqc
//     path '*.pdf'
//
//     script:
//     """
//     RUN_SPP=`which run_spp.R`
//     Rscript -e "library(caTools); source(\\"\$RUN_SPP\\")" -c="${bam[0]}" -savp="${name}.spp.pdf" -savd="${name}.spp.Rdata" -out="${name}.spp.out" -p=$task.cpus
//     cp $spp_correlation_header ${name}_spp_correlation_mqc.tsv
//     Rscript -e "load('${name}.spp.Rdata'); write.table(crosscorr\\\$cross.correlation, file=\\"${name}_spp_correlation_mqc.tsv\\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"
//
//     awk -v OFS='\t' '{print "${name}", \$9}' ${name}.spp.out | cat $spp_nsc_header - > ${name}_spp_nsc_mqc.tsv
//     awk -v OFS='\t' '{print "${name}", \$10}' ${name}.spp.out | cat $spp_rsc_header - > ${name}_spp_rsc_mqc.tsv
//     """
// }
//
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
//  * STEP 6.1: deepTools plotFingerprint
//  */
// process PLOTFINGERPRINT {
//     tag "${ip} vs ${control}"
//     label 'process_high'
//     publishDir "${params.outdir}/bwa/mergedLibrary/deepTools/plotFingerprint", mode: params.publish_dir_mode
//
//     when:
//     !params.skip_plot_fingerprint
//
//     input:
//     tuple val(antibody), val(replicatesExist), val(multipleGroups), val(ip), path(ipbam), val(control), path(controlbam), path(ipflagstat) from ch_group_bam_plotfingerprint
//
//     output:
//     path '*.raw.txt' into ch_plotfingerprint_mqc
//     path '*.{txt,pdf}'
//
//     script:
//     extend = (params.single_end && params.fragment_size > 0) ? "--extendReads ${params.fragment_size}" : ''
//     """
//     plotFingerprint \\
//         --bamfiles ${ipbam[0]} ${controlbam[0]} \\
//         --plotFile ${ip}.plotFingerprint.pdf \\
//         $extend \\
//         --labels $ip $control \\
//         --outRawCounts ${ip}.plotFingerprint.raw.txt \\
//         --outQualityMetrics ${ip}.plotFingerprint.qcmetrics.txt \\
//         --skipZeros \\
//         --JSDsample ${controlbam[0]} \\
//         --numberOfProcessors $task.cpus \\
//         --numberOfSamples $params.fingerprint_bins
//     """
// }
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
