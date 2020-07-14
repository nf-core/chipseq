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
    bwa_base = params.fasta.substring(lastPath+1)
    ch_fasta = file(params.fasta, checkIfExists: true)
} else {
    exit 1, 'Fasta file not specified!'
}

if (params.bwa_index) {
    lastPath = params.bwa_index.lastIndexOf(File.separator)
    bwa_dir  = params.bwa_index.substring(0,lastPath+1)
    bwa_base = params.bwa_index.substring(lastPath+1)
    Channel
        .fromPath(bwa_dir, checkIfExists: true)
        .set { ch_bwa_index }
}

// Save AWS IGenomes file containing annotation version
if (params.anno_readme && file(params.anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(params.anno_readme).copyTo("${params.outdir}/genome/")
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
include { OUTPUT_DOCUMENTATION } from './modules/local/output_documentation' params(params)
include { GET_SOFTWARE_VERSIONS } from './modules/local/get_software_versions' params(params)
include { CHECK_SAMPLESHEET
          get_samplesheet_paths
          get_samplesheet_design } from './modules/local/check_samplesheet' params(params)

// /*
//  * Include nf-core modules
//  */
// include { FASTQC } from './modules/nf-core/fastqc' params(params)
// include { MULTIQC } from './modules/nf-core/multiqc' params(params)

/*
 * Run the workflow
 */
workflow {

    // Get paths to FastQ files
    // [sample, single_end?, [ fastq_1, fastq_2 ]]
    CHECK_SAMPLESHEET(ch_input)
        .reads
        .splitCsv(header:true, sep:',')
        .map { get_samplesheet_paths(it, params.single_end) }
        .set { ch_raw_reads }

    // Get design information from samplesheet
    // [sample, control, antibody, replicatesExist?, multipleGroups?]
    CHECK_SAMPLESHEET
        .out
        .controls
        .splitCsv(header:true, sep:',')
        .map { get_samplesheet_design(it) }
        .set { ch_design }

    // FASTQC(ch_raw_reads_fastqc)

    OUTPUT_DOCUMENTATION(
        ch_output_docs,
        ch_output_docs_images)

    GET_SOFTWARE_VERSIONS()

    // MULTIQC(
    //     summary,
    //     FASTQC.out,
    //     ch_multiqc_config
    // )
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
// /* --                     PREPARE ANNOTATION FILES                        -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
//
// /*
//  * PREPROCESSING: Build BWA index
//  */
// if (!params.bwa_index) {
//     process BWA_INDEX {
//         tag "$fasta"
//         label 'process_high'
//         publishDir path: { params.save_reference ? "${params.outdir}/genome" : params.outdir },
//             saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode
//
//         input:
//         path fasta from ch_fasta
//
//         output:
//         path 'BWAIndex' into ch_bwa_index
//
//         script:
//         """
//         bwa index -a bwtsw $fasta
//         mkdir BWAIndex && mv ${fasta}* BWAIndex
//         """
//     }
// }
//
// /*
//  * PREPROCESSING: Generate gene BED file
//  */
// // If --gtf is supplied along with --genome
// // Make gene bed from supplied --gtf instead of using iGenomes one automatically
// def MAKE_BED = false
// if (!params.gene_bed) {
//     MAKE_BED = true
// } else if (params.genome && params.gtf) {
//     if (params.genomes[ params.genome ].gtf != params.gtf) {
//         MAKE_BED = true
//     }
// }
// if (MAKE_BED) {
//     process MAKE_GENE_BED {
//         tag "$gtf"
//         label 'process_low'
//         publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
//
//         input:
//         path gtf from ch_gtf
//
//         output:
//         path '*.bed' into ch_gene_bed
//
//         script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
//         """
//         gtf2bed $gtf > ${gtf.baseName}.bed
//         """
//     }
// }
//
// /*
//  * PREPROCESSING: Prepare genome intervals for filtering
//  */
// process MAKE_GENOME_FILTER {
//     tag "$fasta"
//     publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
//
//     input:
//     path fasta from ch_fasta
//     path blacklist from ch_blacklist.ifEmpty([])
//
//     output:
//     path "$fasta"                                      // FASTA FILE FOR IGV
//     path '*.fai'                                       // FAI INDEX FOR REFERENCE GENOME
//     path '*.bed' into ch_genome_filter_regions         // BED FILE WITHOUT BLACKLIST REGIONS
//     path '*.sizes' into ch_genome_sizes_bigwig         // CHROMOSOME SIZES FILE FOR BEDTOOLS
//
//     script:
//     blacklist_filter = params.blacklist ? "sortBed -i $blacklist -g ${fasta}.sizes | complementBed -i stdin -g ${fasta}.sizes" : "awk '{print \$1, '0' , \$2}' OFS='\t' ${fasta}.sizes"
//     """
//     samtools faidx $fasta
//     cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
//     $blacklist_filter > ${fasta}.include_regions.bed
//     """
// }
//
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                        FASTQ QC                                     -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
//
// /*
//  * STEP 1: FastQC
//  */
// process FASTQC {
//     tag "$name"
//     label 'process_medium'
//     publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       filename.endsWith('.zip') ? "zips/$filename" : filename
//                 }
//
//     when:
//     !params.skip_fastqc
//
//     input:
//     tuple val(name), path(reads) from ch_raw_reads_fastqc
//
//     output:
//     path '*.{zip,html}' into ch_fastqc_reports_mqc
//
//     script:
//     // Added soft-links to original fastqs for consistent naming in MultiQC
//     if (params.single_end) {
//         """
//         [ ! -f  ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
//         fastqc -q -t $task.cpus ${name}.fastq.gz
//         """
//     } else {
//         """
//         [ ! -f  ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
//         [ ! -f  ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
//         fastqc -q -t $task.cpus ${name}_1.fastq.gz
//         fastqc -q -t $task.cpus ${name}_2.fastq.gz
//         """
//     }
// }
//
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                        ADAPTER TRIMMING                             -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
//
// /*
//  * STEP 2: Trim Galore!
//  */
// if (params.skip_trimming) {
//     ch_trimmed_reads = ch_raw_reads_trimgalore
//     ch_trimgalore_results_mqc = Channel.empty()
//     ch_trimgalore_fastqc_reports_mqc = Channel.empty()
// } else {
//     process TRIMGALORE {
//         tag "$name"
//         label 'process_high'
//         publishDir "${params.outdir}/trim_galore", mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                           if (filename.endsWith('.html')) "fastqc/$filename"
//                           else if (filename.endsWith('.zip')) "fastqc/zips/$filename"
//                           else if (filename.endsWith('trimming_report.txt')) "logs/$filename"
//                           else params.save_trimmed ? filename : null
//                     }
//
//         input:
//         tuple val(name), path(reads) from ch_raw_reads_trimgalore
//
//         output:
//         tuple val(name), path('*.fq.gz') into ch_trimmed_reads
//         path '*.txt' into ch_trimgalore_results_mqc
//         path '*.{zip,html}' into ch_trimgalore_fastqc_reports_mqc
//
//         script:
//         // Calculate number of --cores for TrimGalore based on value of task.cpus
//         // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
//         // See: https://github.com/nf-core/atacseq/pull/65
//         def cores = 1
//         if (task.cpus) {
//             cores = (task.cpus as int) - 4
//             if (params.single_end) cores = (task.cpus as int) - 3
//             if (cores < 1) cores = 1
//             if (cores > 4) cores = 4
//         }
//
//         c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
//         c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
//         tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
//         tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
//         nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''
//
//         // Added soft-links to original fastqs for consistent naming in MultiQC
//         if (params.single_end) {
//             """
//             [ ! -f  ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
//             trim_galore --cores $cores --fastqc --gzip $c_r1 $tpc_r1 $nextseq ${name}.fastq.gz
//             """
//         } else {
//             """
//             [ ! -f  ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
//             [ ! -f  ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
//             trim_galore --cores $cores --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $nextseq ${name}_1.fastq.gz ${name}_2.fastq.gz
//             """
//         }
//     }
// }
//
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                        ALIGN                                        -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
//
// /*
//  * STEP 3.1: Map read(s) with bwa mem
//  */
// process BWA_MEM {
//     tag "$name"
//     label 'process_high'
//
//     input:
//     tuple val(name), path(reads) from ch_trimmed_reads
//     path index from ch_bwa_index.collect()
//
//     output:
//     tuple val(name), path('*.bam') into ch_bwa_bam
//
//     script:
//     prefix = "${name}.Lb"
//     rg = "\'@RG\\tID:${name}\\tSM:${name.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${name}\\tPU:1\'"
//     if (params.seq_center) {
//         rg = "\'@RG\\tID:${name}\\tSM:${name.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${name}\\tPU:1\\tCN:${params.seq_center}\'"
//     }
//     score = params.bwa_min_score ? "-T ${params.bwa_min_score}" : ''
//     """
//     bwa mem \\
//         -t $task.cpus \\
//         -M \\
//         -R $rg \\
//         $score \\
//         ${index}/${bwa_base} \\
//         $reads \\
//         | samtools view -@ $task.cpus -b -h -F 0x0100 -O BAM -o ${prefix}.bam -
//     """
// }
//
// /*
//  * STEP 3.2: Convert BAM to coordinate sorted BAM
//  */
// process SORT_BAM {
//     tag "$name"
//     label 'process_medium'
//     if (params.save_align_intermeds) {
//         publishDir path: "${params.outdir}/bwa/library", mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                           if (filename.endsWith('.flagstat')) "samtools_stats/$filename"
//                           else if (filename.endsWith('.idxstats')) "samtools_stats/$filename"
//                           else if (filename.endsWith('.stats')) "samtools_stats/$filename"
//                           else filename
//                     }
//     }
//
//     input:
//     tuple val(name), path(bam) from ch_bwa_bam
//
//     output:
//     tuple val(name), path('*.sorted.{bam,bam.bai}') into ch_sort_bam_merge
//     path '*.{flagstat,idxstats,stats}' into ch_sort_bam_flagstat_mqc
//
//     script:
//     prefix = "${name}.Lb"
//     """
//     samtools sort -@ $task.cpus -o ${prefix}.sorted.bam -T $name $bam
//     samtools index ${prefix}.sorted.bam
//     samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
//     samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
//     samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
//     """
// }
//
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                    MERGE LIBRARY BAM                                -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
//
// /*
//  * STEP 4.1: Merge BAM files for all libraries from same sample replicate
//  */
// ch_sort_bam_merge
//     .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1] ] }
//     .groupTuple(by: [0])
//     .map { it ->  [ it[0], it[1].flatten() ] }
//     .set { ch_sort_bam_merge }
//
// process MERGED_BAM {
//     tag "$name"
//     label 'process_medium'
//     publishDir "${params.outdir}/bwa/mergedLibrary", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith('.flagstat')) "samtools_stats/$filename"
//                       else if (filename.endsWith('.idxstats')) "samtools_stats/$filename"
//                       else if (filename.endsWith('.stats')) "samtools_stats/$filename"
//                       else if (filename.endsWith('.metrics.txt')) "picard_metrics/$filename"
//                       else params.save_align_intermeds ? filename : null
//                 }
//
//     input:
//     tuple val(name), path(bams) from ch_sort_bam_merge
//
//     output:
//     tuple val(name), path("*${prefix}.sorted.{bam,bam.bai}") into ch_merge_bam_filter,
//                                                                   ch_merge_bam_preseq
//     path '*.{flagstat,idxstats,stats}' into ch_merge_bam_stats_mqc
//     path '*.txt' into ch_merge_bam_metrics_mqc
//
//     script:
//     prefix = "${name}.mLb.mkD"
//     bam_files = bams.findAll { it.toString().endsWith('.bam') }.sort()
//     def avail_mem = 3
//     if (!task.memory) {
//         log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
//     } else {
//         avail_mem = task.memory.toGiga()
//     }
//     if (bam_files.size() > 1) {
//         """
//         picard -Xmx${avail_mem}g MergeSamFiles \\
//             ${'INPUT='+bam_files.join(' INPUT=')} \\
//             OUTPUT=${name}.sorted.bam \\
//             SORT_ORDER=coordinate \\
//             VALIDATION_STRINGENCY=LENIENT \\
//             TMP_DIR=tmp
//         samtools index ${name}.sorted.bam
//
//         picard -Xmx${avail_mem}g MarkDuplicates \\
//             INPUT=${name}.sorted.bam \\
//             OUTPUT=${prefix}.sorted.bam \\
//             ASSUME_SORTED=true \\
//             REMOVE_DUPLICATES=false \\
//             METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt \\
//             VALIDATION_STRINGENCY=LENIENT \\
//             TMP_DIR=tmp
//
//         samtools index ${prefix}.sorted.bam
//         samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
//         samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
//         samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
//         """
//     } else {
//       """
//       picard -Xmx${avail_mem}g MarkDuplicates \\
//           INPUT=${bam_files[0]} \\
//           OUTPUT=${prefix}.sorted.bam \\
//           ASSUME_SORTED=true \\
//           REMOVE_DUPLICATES=false \\
//           METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt \\
//           VALIDATION_STRINGENCY=LENIENT \\
//           TMP_DIR=tmp
//
//       samtools index ${prefix}.sorted.bam
//       samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
//       samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
//       samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
//       """
//     }
// }
//
// /*
//  * STEP 4.2: Filter BAM file at merged library-level
//  */
// process MERGED_BAM_FILTER {
//     tag "$name"
//     label 'process_medium'
//     publishDir path: "${params.outdir}/bwa/mergedLibrary", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (params.single_end || params.save_align_intermeds) {
//                           if (filename.endsWith('.flagstat')) "samtools_stats/$filename"
//                           else if (filename.endsWith('.idxstats')) "samtools_stats/$filename"
//                           else if (filename.endsWith('.stats')) "samtools_stats/$filename"
//                           else if (filename.endsWith('.sorted.bam')) filename
//                           else if (filename.endsWith('.sorted.bam.bai')) filename
//                           else null
//                       }
//                 }
//
//     input:
//     tuple val(name), path(bam) from ch_merge_bam_filter
//     path bed from ch_genome_filter_regions.collect()
//     path bamtools_filter_config from ch_bamtools_filter_config
//
//     output:
//     tuple val(name), path('*.{bam,bam.bai}') into ch_filter_bam
//     tuple val(name), path('*.flagstat') into ch_filter_bam_flagstat
//     path '*.{idxstats,stats}' into ch_filter_bam_stats_mqc
//
//     script:
//     prefix = params.single_end ? "${name}.mLb.clN" : "${name}.mLb.flT"
//     filter_params = params.single_end ? '-F 0x004' : '-F 0x004 -F 0x0008 -f 0x001'
//     dup_params = params.keep_dups ? '' : '-F 0x0400'
//     multimap_params = params.keep_multi_map ? '' : '-q 1'
//     blacklist_params = params.blacklist ? "-L $bed" : ''
//     name_sort_bam = params.single_end ? '' : "samtools sort -n -@ $task.cpus -o ${prefix}.bam -T $prefix ${prefix}.sorted.bam"
//     """
//     samtools view \\
//         $filter_params \\
//         $dup_params \\
//         $multimap_params \\
//         $blacklist_params \\
//         -b ${bam[0]} \\
//         | bamtools filter \\
//             -out ${prefix}.sorted.bam \\
//             -script $bamtools_filter_config
//
//     samtools index ${prefix}.sorted.bam
//     samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
//     samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
//     samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
//
//     $name_sort_bam
//     """
// }
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
