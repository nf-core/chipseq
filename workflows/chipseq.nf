/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    aligners       : [ 'bwa', 'bowtie2', 'star' ]
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowChipseq.initialise(params, log, valid_params)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta,
    params.gtf, params.gff, params.gene_bed,
    params.bwa_index, params.bowtie2_index, params.star_index,
    params.blacklist,
    params.bamtools_filter_pe_config, params.bamtools_filter_se_config
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check alignment parameters
def prepareToolIndices  = []
// if (!params.skip_alignment) { prepareToolIndices << params.aligner        } //DEL

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

// JSON files required by BAMTools for alignment filtering
ch_bamtools_filter_se_config = file(params.bamtools_filter_se_config, checkIfExists: true)
ch_bamtools_filter_pe_config = file(params.bamtools_filter_pe_config, checkIfExists: true)

// Header files for MultiQC
ch_spp_nsc_header           = file("$projectDir/assets/multiqc/spp_nsc_header.txt", checkIfExists: true)
ch_spp_rsc_header           = file("$projectDir/assets/multiqc/spp_rsc_header.txt", checkIfExists: true)
ch_spp_correlation_header   = file("$projectDir/assets/multiqc/spp_correlation_header.txt", checkIfExists: true)
ch_peak_count_header        = file("$projectDir/assets/multiqc/peak_count_header.txt", checkIfExists: true)
ch_frip_score_header        = file("$projectDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)
ch_peak_annotation_header   = file("$projectDir/assets/multiqc/peak_annotation_header.txt", checkIfExists: true)
ch_deseq2_pca_header        = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_deseq2_clustering_header = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BEDTOOLS_GENOMECOV                  } from '../modules/local/bedtools_genomecov'
include { FRIP_SCORE                          } from '../modules/local/frip_score'
include { PLOT_MACS2_QC                       } from '../modules/local/plot_macs2_qc'
include { PLOT_HOMER_ANNOTATEPEAKS            } from '../modules/local/plot_homer_annotatepeaks'
include { MACS2_CONSENSUS                     } from '../modules/local/macs2_consensus'
include { DESEQ2_QC                           } from '../modules/local/deseq2_qc'
include { IGV                                 } from '../modules/local/igv'
include { MULTIQC                             } from '../modules/local/multiqc'
include { MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS } from '../modules/local/multiqc_custom_phantompeakqualtools'
include { MULTIQC_CUSTOM_PEAKS                } from '../modules/local/multiqc_custom_peaks'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK         } from '../subworkflows/local/input_check'
include { PREPARE_GENOME      } from '../subworkflows/local/prepare_genome'
include { FILTER_BAM_BAMTOOLS } from '../subworkflows/local/filter_bam_bamtools'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { PICARD_MERGESAMFILES          } from '../modules/nf-core/modules/picard/mergesamfiles/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from '../modules/nf-core/modules/picard/collectmultiplemetrics/main'
include { PRESEQ_LCEXTRAP               } from '../modules/nf-core/modules/preseq/lcextrap/main'
include { PHANTOMPEAKQUALTOOLS          } from '../modules/nf-core/modules/phantompeakqualtools/main'
include { UCSC_BEDGRAPHTOBIGWIG         } from '../modules/nf-core/modules/ucsc/bedgraphtobigwig/main'
include { DEEPTOOLS_COMPUTEMATRIX       } from '../modules/nf-core/modules/deeptools/computematrix/main'
include { DEEPTOOLS_PLOTPROFILE         } from '../modules/nf-core/modules/deeptools/plotprofile/main'
include { DEEPTOOLS_PLOTHEATMAP         } from '../modules/nf-core/modules/deeptools/plotheatmap/main'
include { DEEPTOOLS_PLOTFINGERPRINT     } from '../modules/nf-core/modules/deeptools/plotfingerprint/main'
include { MACS2_CALLPEAK                } from '../modules/nf-core/modules/macs2/callpeak/main'
include { SUBREAD_FEATURECOUNTS         } from '../modules/nf-core/modules/subread/featurecounts/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_MACS2     } from '../modules/nf-core/modules/homer/annotatepeaks/main'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_CONSENSUS } from '../modules/nf-core/modules/homer/annotatepeaks/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//

include { FASTQC_TRIMGALORE      } from '../subworkflows/nf-core/fastqc_trimgalore'
include { ALIGN_BWA_MEM          } from '../subworkflows/nf-core/align_bwa_mem'
include { ALIGN_BOWTIE2          } from '../subworkflows/nf-core/align_bowtie2'
include { ALIGN_STAR             } from '../subworkflows/nf-core/align_star'
include { MARK_DUPLICATES_PICARD } from '../subworkflows/nf-core/mark_duplicates_picard'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CHIPSEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME (
        params.aligner
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input),
        params.seq_center
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        INPUT_CHECK.out.reads,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

    //
    // SUBWORKFLOW: Alignment with BWA & BAM QC
    //
    ch_genome_bam        = Channel.empty()
    ch_genome_bam_index  = Channel.empty()
    ch_samtools_stats    = Channel.empty()
    ch_samtools_flagstat = Channel.empty()
    ch_samtools_idxstats = Channel.empty()
    if (params.aligner == 'bwa') {
        ALIGN_BWA_MEM (
            FASTQC_TRIMGALORE.out.reads,
            PREPARE_GENOME.out.bwa_index
        )
        ch_genome_bam        = ALIGN_BWA_MEM.out.bam
        ch_genome_bam_index  = ALIGN_BWA_MEM.out.bai
        ch_samtools_stats    = ALIGN_BWA_MEM.out.stats
        ch_samtools_flagstat = ALIGN_BWA_MEM.out.flagstat
        ch_samtools_idxstats = ALIGN_BWA_MEM.out.idxstats
        ch_versions = ch_versions.mix(ALIGN_BWA_MEM.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with BOWTIE2 & BAM QC
    //
    if (params.aligner == 'bowtie2') {
        ALIGN_BOWTIE2 (
            FASTQC_TRIMGALORE.out.reads,
            PREPARE_GENOME.out.bowtie2_index,
            params.save_unaligned
        )
        ch_genome_bam        = ALIGN_BOWTIE2.out.bam
        ch_genome_bam_index  = ALIGN_BOWTIE2.out.bai
        ch_samtools_stats    = ALIGN_BOWTIE2.out.stats
        ch_samtools_flagstat = ALIGN_BOWTIE2.out.flagstat
        ch_samtools_idxstats = ALIGN_BOWTIE2.out.idxstats
        ch_versions = ch_versions.mix(ALIGN_BOWTIE2.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with STAR & BAM QC
    //
    if (params.aligner == 'star') {
        ALIGN_STAR (
            FASTQC_TRIMGALORE.out.reads,
            PREPARE_GENOME.out.star_index
        )
        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bam_index  = ALIGN_STAR.out.bai
        ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final

        ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)
    }

    //
    // SUBWORKFLOW: Merge resequenced BAM files
    //
    ch_genome_bam
        .map {
            meta, bam ->
                fmeta = meta.findAll { it.key != 'read_group' }
                fmeta.id = fmeta.id.split('_')[0..-2].join('_')
                [ fmeta, bam ] }
        .groupTuple(by: [0])
        .map { it ->  [ it[0], it[1].flatten() ] }
        .set { ch_sort_bam }

    PICARD_MERGESAMFILES (
        ch_sort_bam
    )
    ch_versions = ch_versions.mix(PICARD_MERGESAMFILES.out.versions.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Mark duplicates & filter BAM files after merging
    //
    ch_markduplicates_multiqc = Channel.empty()
    MARK_DUPLICATES_PICARD (
        PICARD_MERGESAMFILES.out.bam
    )

    //
    // SUBWORKFLOW: Fix getting name sorted BAM here for PE/SE
    //
    FILTER_BAM_BAMTOOLS (
        MARK_DUPLICATES_PICARD.out.bam.join(MARK_DUPLICATES_PICARD.out.bai, by: [0]),
        PREPARE_GENOME.out.filtered_bed,

        ch_bamtools_filter_se_config,
        ch_bamtools_filter_pe_config
    )
    ch_versions = ch_versions.mix(FILTER_BAM_BAMTOOLS.out.versions.first().ifEmpty(null))

    //
    // MODULE: Post alignment QC
    //
    ch_picardcollectmultiplemetrics_multiqc = Channel.empty()
    if (!params.skip_picard_metrics) {
        PICARD_COLLECTMULTIPLEMETRICS (
            FILTER_BAM_BAMTOOLS.out.bam,
            PREPARE_GENOME.out.fasta
        )
        ch_picardcollectmultiplemetrics_multiqc = MARK_DUPLICATES_PICARD.out.metrics
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
    }

    //
    // MODULE: Library coverage
    //
    ch_preseq_multiqc = Channel.empty()
    if (!params.skip_preseq) {
        PRESEQ_LCEXTRAP (
            FILTER_BAM_BAMTOOLS.out.bam
        )
        ch_preseq_multiqc = PRESEQ_LCEXTRAP.out.ccurve
        ch_versions       = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())
    }

    //
    // MODULE: Strand cross-correlation
    //
    PHANTOMPEAKQUALTOOLS (
        FILTER_BAM_BAMTOOLS.out.bam
    )
    ch_versions = ch_versions.mix(PHANTOMPEAKQUALTOOLS.out.versions.first())

    MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS (
        PHANTOMPEAKQUALTOOLS.out.spp.join(PHANTOMPEAKQUALTOOLS.out.rdata, by: [0]),
        ch_spp_nsc_header,
        ch_spp_rsc_header,
        ch_spp_correlation_header
    )

    //
    // MODULE: Coverage tracks
    //
    BEDTOOLS_GENOMECOV (
        FILTER_BAM_BAMTOOLS.out.bam.join(FILTER_BAM_BAMTOOLS.out.flagstat, by: [0])
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

    //
    // MODULE: Coverage tracks
    //
    UCSC_BEDGRAPHTOBIGWIG (
        BEDTOOLS_GENOMECOV.out.bedgraph,
        PREPARE_GENOME.out.chrom_sizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())

    //
    // MODULE: Coverage plots
    //
    ch_deeptoolsplotprofile_multiqc = Channel.empty()
    if (!params.skip_plot_profile) {
        DEEPTOOLS_COMPUTEMATRIX (
            UCSC_BEDGRAPHTOBIGWIG.out.bigwig,
            PREPARE_GENOME.out.gene_bed
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX.out.versions.first())

        DEEPTOOLS_PLOTPROFILE (
            DEEPTOOLS_COMPUTEMATRIX.out.matrix
        )
        ch_deeptoolsplotprofile_multiqc = DEEPTOOLS_PLOTPROFILE.out.table
        ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX.out.versions.first())

        DEEPTOOLS_PLOTHEATMAP (
            DEEPTOOLS_COMPUTEMATRIX.out.matrix
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTHEATMAP.out.versions.first())
    }

    //
    // Refactor channels: [ val(meta), [ ip_bam, control_bam ] [ ip_bai, control_bai ] ]
    //
    FILTER_BAM_BAMTOOLS
        .out
        .bam
        .join (FILTER_BAM_BAMTOOLS.out.bai, by: [0])
        .map {
            meta, bam, bai ->
                meta.control ? null : [ meta.id, [ bam ] , [ bai ] ]
        }
        .set { ch_control_bam_bai }

    FILTER_BAM_BAMTOOLS
        .out
        .bam
        .join (FILTER_BAM_BAMTOOLS.out.bai, by: [0])
        .map {
            meta, bam, bai ->
                meta.control ? [ meta.control, meta, [ bam ], [ bai ] ] : null
        }
        .combine(ch_control_bam_bai, by: 0)
        .map { it -> [ it[1] , it[2] + it[4], it[3] + it[5] ] }
        .set { ch_ip_control_bam_bai }

    //
    // plotFingerprint for IP and control together
    //
    ch_deeptoolsplotfingerprint_multiqc = Channel.empty()
    if (!params.skip_plot_fingerprint) {
        DEEPTOOLS_PLOTFINGERPRINT (
            ch_ip_control_bam_bai
        )
        ch_deeptools_plotfingerprintmultiqc = DEEPTOOLS_PLOTFINGERPRINT.out.matrix
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTFINGERPRINT.out.versions.first())
    }

    //
    // Call peaks
    //
    ch_custompeaks_frip_multiqc       = Channel.empty()
    ch_custompeaks_count_multiqc      = Channel.empty()
    ch_plothomerannotatepeaks_multiqc = Channel.empty()
    ch_subreadfeaturecounts_multiqc   = Channel.empty()
    if (params.macs_gsize) {
        // Create channel: [ val(meta), ip_bam, control_bam ]
        ch_ip_control_bam_bai
            .map { meta, bams, bais -> [ meta , bams[0], bams[1] ] }
            .set { ch_ip_control_bam }


        MACS2_CALLPEAK (
            ch_ip_control_bam,
            params.macs_gsize
        )
        ch_versions = ch_versions.mix(MACS2_CALLPEAK.out.versions.first())

        ch_ip_control_bam
            .join(MACS2_CALLPEAK.out.peak, by: [0])
            .map { it -> [ it[0], it[1], it[3] ] }
            .set { ch_ip_peak }
        FRIP_SCORE (
            ch_ip_peak
        )
        ch_versions = ch_versions.mix(FRIP_SCORE.out.versions.first())

        ch_ip_peak
            .join(FRIP_SCORE.out.txt, by: [0])
            .map { it -> [ it[0], it[2], it[3] ] }
            .set { ch_ip_peak_frip }

        MULTIQC_CUSTOM_PEAKS (
            ch_ip_peak_frip,
            ch_peak_count_header,
            ch_frip_score_header
        )
        ch_custompeaks_frip_multiqc  = MULTIQC_CUSTOM_PEAKS.out.frip
        ch_custompeaks_count_multiqc = MULTIQC_CUSTOM_PEAKS.out.count

        if (!params.skip_peak_annotation && !params.skip_peak_qc) {
            PLOT_MACS2_QC (
                MACS2_CALLPEAK.out.peak.collect{it[1]}
            )
            ch_versions = ch_versions.mix(PLOT_MACS2_QC.out.versions)

            HOMER_ANNOTATEPEAKS_MACS2 (
                MACS2_CALLPEAK.out.peak,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.gtf
            )
            ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS_MACS2.out.versions.first())

            PLOT_HOMER_ANNOTATEPEAKS (
                HOMER_ANNOTATEPEAKS_MACS2.out.txt.collect{it[1]},
                ch_peak_annotation_header,
                "_peaks.annotatePeaks.txt"
            )
            ch_plothomerannotatepeaks_multiqc = PLOT_HOMER_ANNOTATEPEAKS.out.tsv
            ch_versions = ch_versions.mix(PLOT_HOMER_ANNOTATEPEAKS.out.versions)
        }

        //
        //  Consensus peaks analysis
        //
        if (!params.skip_consensus_peaks) {
            // Create channel: [ meta , [ peaks ] ]
            // Where meta = [ id:antibody, multiple_groups:true/false, replicates_exist:true/false ]
            MACS2_CALLPEAK
                .out
                .peak
                .map { meta, peak -> [ meta.antibody, meta.id.split('_')[0..-2].join('_'), peak ] }
                .groupTuple()
                .map {
                    antibody, groups, peaks ->
                        [
                            antibody,
                            groups.groupBy().collectEntries { [(it.key) : it.value.size()] },
                            peaks
                        ] }
                .map {
                    antibody, groups, peaks ->
                        def meta = [:]
                        meta.id = antibody
                        meta.multiple_groups = groups.size() > 1
                        meta.replicates_exist = groups.max { groups.value }.value > 1
                        [ meta, peaks ] }
                .set { ch_antibody_peaks }

            MACS2_CONSENSUS (
                ch_antibody_peaks
            )
            ch_versions = ch_versions.mix(MACS2_CONSENSUS.out.versions)

            if (!params.skip_peak_annotation) {
                HOMER_ANNOTATEPEAKS_CONSENSUS (
                    MACS2_CONSENSUS.out.bed,
                    PREPARE_GENOME.out.fasta,
                    PREPARE_GENOME.out.gtf
                )
                ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS_CONSENSUS.out.versions)
                // cut -f2- ${prefix}.annotatePeaks.txt | awk 'NR==1; NR > 1 {print \$0 | "sort -T '.' -k1,1 -k2,2n"}' | cut -f6- > tmp.txt
                // paste $bool tmp.txt > ${prefix}.boolean.annotatePeaks.txt
            }

            // Create channel: [ val(meta), ip_bam ]
            MACS2_CONSENSUS
                .out
                .saf
                .map { meta, saf -> [ meta.id, meta, saf ] }
                .set { ch_ip_saf }

            ch_ip_control_bam
                .map { meta, ip_bam, control_bam -> [ meta.antibody, meta, ip_bam ] }
                .groupTuple()
                .map { it -> [ it[0], it[1][0], it[2].flatten().sort() ] }
                .join(ch_ip_saf)
                .map {
                    it ->
                        fmeta = it[1]
                        fmeta['id'] = it[3]['id']
                        fmeta['replicates_exist'] = it[3]['replicates_exist']
                        fmeta['multiple_groups']  = it[3]['multiple_groups']
                        [ fmeta, it[2], it[4] ] }
                .set { ch_ip_bam }

            SUBREAD_FEATURECOUNTS (
                ch_ip_bam
            )
            ch_subreadfeaturecounts_multiqc = SUBREAD_FEATURECOUNTS.out.summary
            ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

            ch_deseq2_pca_header = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
            ch_deseq2_clustering_header = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
            if (!params.skip_deseq2_qc) {
                DESEQ2_QC (
                    SUBREAD_FEATURECOUNTS.out.counts,
                    ch_deseq2_pca_header,
                    ch_deseq2_clustering_header
                )
            }

            //
            // Create IGV session
            //
            if (!params.skip_igv) {
                IGV (
                    PREPARE_GENOME.out.fasta,
                    UCSC_BEDGRAPHTOBIGWIG.out.bigwig.collect{it[1]}.ifEmpty([]),
                    MACS2_CALLPEAK.out.peak.collect{it[1]}.ifEmpty([]),
                    MACS2_CONSENSUS.out.bed.collect{it[1]}.ifEmpty([]),
                    "bwa/mergedLibrary/bigwig",
                    { ["bwa/mergedLibrary/macs2",
                        params.narrow_peak? '/narrowPeak' : '/broadPeak'
                        ].join('') },
                    { ["bwa/mergedLibrary/macs2",
                        params.narrow_peak? '/narrowPeak' : '/broadPeak'
                        ].join('') }
                )
                ch_versions = ch_versions.mix(IGV.out.versions)
            }
        }
    }

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowChipseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),

            FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),

            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),

            MARK_DUPLICATES_PICARD.out.stats.collect{it[1]}.ifEmpty([]),
            MARK_DUPLICATES_PICARD.out.flagstat.collect{it[1]}.ifEmpty([]),
            MARK_DUPLICATES_PICARD.out.idxstats.collect{it[1]}.ifEmpty([]),
            MARK_DUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]),

            FILTER_BAM_BAMTOOLS.out.stats.collect{it[1]}.ifEmpty([]),
            FILTER_BAM_BAMTOOLS.out.flagstat.collect{it[1]}.ifEmpty([]),
            FILTER_BAM_BAMTOOLS.out.idxstats.collect{it[1]}.ifEmpty([]),
            ch_picardcollectmultiplemetrics_multiqc.collect{it[1]}.ifEmpty([]),

            ch_preseq_multiqc.collect{it[1]}.ifEmpty([]),
            ch_deeptoolsplotprofile_multiqc.collect{it[1]}.ifEmpty([]),
            ch_deeptoolsplotfingerprint_multiqc.collect{it[1]}.ifEmpty([]),
            PHANTOMPEAKQUALTOOLS.out.spp.collect{it[1]}.ifEmpty([]),
            MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.nsc.collect{it[1]}.ifEmpty([]),
            MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.rsc.collect{it[1]}.ifEmpty([]),
            MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.correlation.collect{it[1]}.ifEmpty([]),

            ch_custompeaks_frip_multiqc.collect{it[1]}.ifEmpty([]),
            ch_custompeaks_count_multiqc.collect{it[1]}.ifEmpty([]),
            ch_plothomerannotatepeaks_multiqc.collect{it[1]}.ifEmpty([]),
            ch_subreadfeaturecounts_multiqc.collect{it[1]}.ifEmpty([])//,
            // path ('macs/consensus/*') from ch_macs_consensus_deseq_mqc.collect().ifEmpty([])
        )
        multiqc_report       = MULTIQC.out.report.toList()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
