/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowChipseq.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta,
    params.gtf, params.gff, params.gene_bed,
    params.bwa_index,
    params.blacklist,
    params.bamtools_filter_pe_config, params.bamtools_filter_se_config
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check blacklist file
ch_blacklist = params.blacklist ? file(params.blacklist) : Channel.empty()

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
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

def peakType = params.narrow_peak ? 'narrowPeak' : 'broadPeak'

/*
========================================================================================
    CONFIG FILES
========================================================================================
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
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def plot_homer_annotatepeaks_options          = modules['plot_homer_annotatepeaks']
plot_homer_annotatepeaks_options.publish_dir  += "/$peakType/qc"

def plot_macs2_qc_options         = modules['plot_macs2_qc']
plot_macs2_qc_options.publish_dir += "/$peakType/qc"

def macs2_consensus_options         = modules['macs2_consensus']
macs2_consensus_options.publish_dir += "/$peakType/qc"

def multiqc_options         = modules['multiqc']
multiqc_options.args        += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''
multiqc_options.publish_dir += "/$peakType"

include { GTF2BED                             } from '../modules/local/gtf2bed'                              addParams( options: [:] )
include { GET_CHROM_SIZES                     } from '../modules/local/get_chrom_sizes'                      addParams( options: [:] )
include { MAKE_GENOME_FILTER                  } from '../modules/local/make_genome_filter'                   addParams( options: [:] )
include { BEDTOOLS_GENOMECOV                  } from '../modules/local/bedtools_genomecov'                   addParams( options: modules['bedtools_genomecov'] )
include { PLOT_HOMER_ANNOTATEPEAKS            } from '../modules/local/plot_homer_annotatepeaks'             addParams( options: plot_homer_annotatepeaks_options )
include { PLOT_MACS2_QC                       } from '../modules/local/plot_macs2_qc'                        addParams( options: plot_macs2_qc_options )
include { MACS2_CONSENSUS                     } from '../modules/local//macs2_consensus'                     addParams( options: macs2_consensus_options )
include { FRIP_SCORE                          } from '../modules/local/frip_score'                           addParams( options: modules['frip_score'] )
//include { DESEQ2_QC  } from '../modules/local/deseq2_qc'                             addParams( options: deseq2_qc_options, multiqc_label: 'star_salmon'   )
include { IGV                                 } from '../modules/local/igv'                                  addParams( options: [:] )
include { OUTPUT_DOCUMENTATION                } from '../modules/local/output_documentation'                 addParams( options: [:] )
include { MULTIQC                             } from '../modules/local/multiqc'                              addParams( options: multiqc_options )
include { MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS } from '../modules/local/multiqc_custom_phantompeakqualtools'  addParams( options: modules['multiqc_custom_phantompeakqualtools'] )
include { MULTIQC_CUSTOM_PEAKS                } from '../modules/local/multiqc_custom_peaks'                 addParams( options: modules['multiqc_custom_peaks'] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )
include { BAM_CLEAN   } from '../subworkflows/local/bam_clean'   addParams( bam_filter_options: modules['bam_filter'], bam_remove_orphans_options: modules['bam_remove_orphans'], samtools_sort_options: modules['samtools_sort_filter'], samtools_index_options:  modules['samtools_sort_filter'], samtools_stats_options:  modules['samtools_sort_filter'] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
def deeptools_plotfingerprint_options = modules['deeptools_plotfingerprint']
deeptools_plotfingerprint_options.args += " --numberOfSamples $params.fingerprint_bins"

def macs2_callpeak_options            = modules['macs2_callpeak']
macs2_callpeak_options.args           += params.narrow_peak ? '' : Utils.joinModuleArgs(['--broad', "--broad-cutoff ${params.broad_cutoff}"])
macs2_callpeak_options.args           += params.save_macs_pileup ? Utils.joinModuleArgs(['--bdg', '--SPMR']) : ''
macs2_callpeak_options.args           += params.macs_fdr ? Utils.joinModuleArgs(["--qvalue ${params.macs_fdr}"]) : ''
macs2_callpeak_options.args           += params.macs_pvalue ? Utils.joinModuleArgs(["--pvalue ${params.macs_pvalue}"]) : ''
macs2_callpeak_options['publish_dir'] += "/$peakType"

def subread_featurecounts_options            = modules['subread_featurecounts']
subread_featurecounts_options['publish_dir'] += "/$peakType/consensus"

def homer_annotatepeaks_macs2_options            = modules['homer_annotatepeaks_macs2']
homer_annotatepeaks_macs2_options['publish_dir'] += "/$peakType"

def homer_annotatepeaks_consensus_options            = modules['homer_annotatepeaks_consensus']
homer_annotatepeaks_consensus_options['publish_dir'] += "/$peakType/consensus"

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
// include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   ) // TODO review
include { BWA_INDEX                     } from '../modules/nf-core/modules/bwa/index/main'                     addParams( options: modules['bwa_index'] )
include { PICARD_MERGESAMFILES          } from '../modules/nf-core/modules/picard/mergesamfiles/main'          addParams( options: modules['picard_mergesamfiles'] )
include { PICARD_COLLECTMULTIPLEMETRICS } from '../modules/nf-core/modules/picard/collectmultiplemetrics/main' addParams( options: modules['picard_collectmultiplemetrics'] )
include { PRESEQ_LCEXTRAP               } from '../modules/nf-core/modules/preseq/lcextrap/main'               addParams( options: modules['preseq_lcextrap'] )
include { UCSC_BEDGRAPHTOBIGWIG         } from '../modules/nf-core/modules/ucsc/bedgraphtobigwig/main'         addParams( options: modules['ucsc_bedgraphtobigwig'] )
include { DEEPTOOLS_COMPUTEMATRIX       } from '../modules/nf-core/modules/deeptools/computematrix/main'       addParams( options: modules['deeptools_computematrix'] )
include { DEEPTOOLS_PLOTPROFILE         } from '../modules/nf-core/modules/deeptools/plotprofile/main'         addParams( options: modules['deeptools_plotprofile'] )
include { DEEPTOOLS_PLOTHEATMAP         } from '../modules/nf-core/modules/deeptools/plotheatmap/main'         addParams( options: modules['deeptools_plotheatmap'] )
include { DEEPTOOLS_PLOTFINGERPRINT     } from '../modules/nf-core/modules/deeptools/plotfingerprint/main'     addParams( options: deeptools_plotfingerprint_options )
include { PHANTOMPEAKQUALTOOLS          } from '../modules/nf-core/modules/phantompeakqualtools/main'          addParams( options: modules['phantompeakqualtools'] )
include { MACS2_CALLPEAK                } from '../modules/nf-core/modules/macs2/callpeak/main'                addParams( options: macs2_callpeak_options )
include { SUBREAD_FEATURECOUNTS         } from '../modules/nf-core/modules/subread/featurecounts/main'         addParams( options: subread_featurecounts_options )
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'   addParams( options: [publish_files : ['_versions.yml':'']] )

include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_MACS2 }     from '../modules/nf-core/modules/homer/annotatepeaks/main' addParams( options:homer_annotatepeaks_macs2_options )
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_CONSENSUS } from '../modules/nf-core/modules/homer/annotatepeaks/main' addParams( options:homer_annotatepeaks_consensus_options )

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { FASTQC_TRIMGALORE             } from '../subworkflows/nf-core/fastqc_trimgalore'      addParams( fastqc_options: modules['fastqc'], trimgalore_options: modules['trimgalore'] )
include { MAP_BWA_MEM                   } from '../subworkflows/nf-core/map_bwa_mem'            addParams( bwa_mem_options: modules['bwa_mem'], samtools_sort_options: modules['samtools_sort_lib'], samtools_index_options:  modules['samtools_sort_lib'], samtools_stats_options:  modules['samtools_sort_lib'] )
include { MARK_DUPLICATES_PICARD        } from '../subworkflows/nf-core/mark_duplicates_picard' addParams( markduplicates_options: modules['picard_markduplicates'], samtools_index_options: modules['samtools_sort_merged_lib'] , samtools_stats_options: modules['samtools_sort_merged_lib']  )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CHIPSEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input),
        params.seq_center,
        [:]
    )

    //
    // MODULE: Run FastQC
    //
    // TODO check if this fastqc process from nf-core is needed Don't think so done in FASTQC_TRIMGALORE
    // FASTQC (
    //     INPUT_CHECK.out.reads
    // )
    // ch_software_versions = ch_software_versions.mix(FASTQC.out.versions.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Prepare genome files
    //
    ch_index = params.bwa_index ? Channel.value(file(params.bwa_index)) : BWA_INDEX ( ch_fasta ).index

    if (makeBED) { ch_gene_bed = GTF2BED ( ch_gtf ) }

    MAKE_GENOME_FILTER (
        GET_CHROM_SIZES ( ch_fasta ).sizes,
        ch_blacklist.ifEmpty([])
    )
    // ch_software_versions = Channel.empty()
    // ch_software_versions = ch_software_versions.mix(MAKE_GENOME_FILTER.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(MAKE_GENOME_FILTER.out.versions)

    //
    // SUBWORKFLOW: Read QC & trimming
    //
    nextseq = params.trim_nextseq > 0 ? " --nextseq ${params.trim_nextseq}" : ''
    params.modules['trimgalore'].args += nextseq
    FASTQC_TRIMGALORE (
        INPUT_CHECK.out.reads,
        params.skip_fastqc,
        params.skip_trimming
    )
    // ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.fastqc_versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

    //
    // SUBWORKFLOW: Map reads & BAM QC
    //
    score = params.bwa_min_score ? " -T ${params.bwa_min_score}" : ''
    params.modules['bwa_mem'].args += score
    MAP_BWA_MEM (
        FASTQC_TRIMGALORE.out.reads,
        ch_index
    )
    // ch_software_versions = ch_software_versions.mix(MAP_BWA_MEM.out.bwa_versions.first())
    // ch_software_versions = ch_software_versions.mix(MAP_BWA_MEM.out.samtools_versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(MAP_BWA_MEM.out.versions)

    //
    // SUBWORKFLOW: Merge resequenced BAM files
    //
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
        ch_sort_bam
    )
    // ch_software_versions = ch_software_versions.mix(PICARD_MERGESAMFILES.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(PICARD_MERGESAMFILES.out.versions)

    //
    // SUBWORKFLOW: Mark duplicates & filter BAM files
    //
    MARK_DUPLICATES_PICARD (
        PICARD_MERGESAMFILES.out.bam
    )

    //
    // SUBWORKFLOW: Fix getting name sorted BAM here for PE/SE
    //
    BAM_CLEAN (
        MARK_DUPLICATES_PICARD.out.bam.join(MARK_DUPLICATES_PICARD.out.bai, by: [0]),
        MAKE_GENOME_FILTER.out.bed.collect(),
        ch_bamtools_filter_se_config,
        ch_bamtools_filter_pe_config
    )
    // ch_software_versions = ch_software_versions.mix(BAM_CLEAN.out.bamtools_versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(BAM_CLEAN.out.versions)

    //
    // MODULE: Post alignment QC
    //
    PICARD_COLLECTMULTIPLEMETRICS (
        BAM_CLEAN.out.bam,
        ch_fasta
    )

    // TODO Add versions to picard collectmultiplemetrics

    //
    // MODULE: Library coverage
    //
    PRESEQ_LCEXTRAP (
        BAM_CLEAN.out.bam
    )
    // ch_software_versions = ch_software_versions.mix(PRESEQ_LCEXTRAP.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions)

    //
    // MODULE: Strand cross-correlation
    //
    PHANTOMPEAKQUALTOOLS (
        BAM_CLEAN.out.bam
    )
    // ch_software_versions = ch_software_versions.mix(PHANTOMPEAKQUALTOOLS.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(PHANTOMPEAKQUALTOOLS.out.versions)

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
        BAM_CLEAN.out.bam.join(BAM_CLEAN.out.flagstat, by: [0])
    )

    //
    // MODULE: Coverage tracks
    //
    UCSC_BEDGRAPHTOBIGWIG (
        BEDTOOLS_GENOMECOV.out.bedgraph,
        GET_CHROM_SIZES.out.sizes
    )
    // ch_software_versions = ch_software_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions)

    //
    // MODULE: Coverage plots
    //
    DEEPTOOLS_COMPUTEMATRIX (
        UCSC_BEDGRAPHTOBIGWIG.out.bigwig,
        ch_gene_bed
    )
    // ch_software_versions = ch_software_versions.mix(DEEPTOOLS_COMPUTEMATRIX.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX.out.versions)

    DEEPTOOLS_PLOTPROFILE (
        DEEPTOOLS_COMPUTEMATRIX.out.matrix
    )

    DEEPTOOLS_PLOTHEATMAP (
        DEEPTOOLS_COMPUTEMATRIX.out.matrix
    )

    //
    // Refactor channels: [ val(meta), [ ip_bam, control_bam ] [ ip_bai, control_bai ] ]
    //
    BAM_CLEAN
        .out
        .bam
        .join ( BAM_CLEAN.out.bai, by: [0] )
        .map { meta, bam, bai -> meta.control ? null : [ meta.id, [ bam ] , [ bai ] ] }
        .set { ch_control_bam_bai }

    BAM_CLEAN
        .out
        .bam
        .join ( BAM_CLEAN.out.bai, by: [0] )
        .map { meta, bam, bai -> meta.control ? [ meta.control, meta, [ bam ], [ bai ] ] : null }
        .combine(ch_control_bam_bai, by: 0)
        .map { it -> [ it[1] , it[2] + it[4], it[3] + it[5] ] }
        .set { ch_ip_control_bam_bai }

    //
    // plotFingerprint for IP and control together
    //
    DEEPTOOLS_PLOTFINGERPRINT (
        ch_ip_control_bam_bai
    )

    if (params.macs_gsize) {

        //
        // Call peaks
        //

        // Create channel: [ val(meta), ip_bam, control_bam ]
        ch_ip_control_bam_bai
            .map { meta, bams, bais -> [ meta , bams[0], bams[1] ] }
            .set { ch_ip_control_bam }

        MACS2_CALLPEAK (
            ch_ip_control_bam,
            params.macs_gsize
        )
        // ch_software_versions = ch_software_versions.mix(MACS2_CALLPEAK.out.versions.first().ifEmpty(null))
        ch_versions = ch_versions.mix(MACS2_CALLPEAK.out.versions)

        ch_ip_control_bam
            .join(MACS2_CALLPEAK.out.peak, by: [0])
            .map { it -> [ it[0], it[1], it[3] ] }
            .set { ch_ip_peak }
        FRIP_SCORE (
            ch_ip_peak
        )

        ch_ip_peak
            .join(FRIP_SCORE.out.txt, by: [0])
            .map { it -> [ it[0], it[2], it[3] ] }
            .set { ch_ip_peak_frip }
        MULTIQC_CUSTOM_PEAKS (
            ch_ip_peak_frip,
            ch_peak_count_header,
            ch_frip_score_header
        )

        PLOT_MACS2_QC (
            MACS2_CALLPEAK.out.peak.collect{it[1]}
        )

        HOMER_ANNOTATEPEAKS_MACS2 (
            MACS2_CALLPEAK.out.peak,
            ch_fasta,
            ch_gtf
        )
        // ch_software_versions = ch_software_versions.mix(HOMER_ANNOTATEPEAKS_MACS2.out.versions.first().ifEmpty(null))
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS_MACS2.out.versions)

        PLOT_HOMER_ANNOTATEPEAKS (
            HOMER_ANNOTATEPEAKS_MACS2.out.txt.collect{it[1]},
            ch_peak_annotation_header,
            "_peaks.annotatePeaks.txt"
        )

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

        HOMER_ANNOTATEPEAKS_CONSENSUS (
            MACS2_CONSENSUS.out.bed,
            ch_fasta,
            ch_gtf
        )
        // cut -f2- ${prefix}.annotatePeaks.txt | awk 'NR==1; NR > 1 {print \$0 | "sort -T '.' -k1,1 -k2,2n"}' | cut -f6- > tmp.txt
        // paste $bool tmp.txt > ${prefix}.boolean.annotatePeaks.txt

        // Create channel: [ val(meta), ip_bam ]
        MACS2_CONSENSUS
            .out
            .saf
            .map { meta, saf -> [ meta.id, meta, saf ] }
            .set { ch_ip_saf }

        ch_ip_control_bam
            .map { meta, ip_bam, control_bam -> [ meta.antibody, meta, ip_bam ] }
            .combine(ch_ip_saf)
            .map {
                it ->
                    fmeta = it[1]
                    fmeta['replicates_exist'] = it[4]['replicates_exist']
                    fmeta['multiple_groups'] = it[4]['multiple_groups']
                    [ fmeta, it[2], it[5] ] }
            .set { ch_ip_bam }

        SUBREAD_FEATURECOUNTS (
            ch_ip_bam
        )
        // ch_software_versions = ch_software_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first().ifEmpty(null))
        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions)

        //
        // DESEQ2_FEATURECOUNTS (
        //     params.modules['deseq2_featurecounts']
        // )
        // ch_deseq2_pca_header = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
        // ch_deseq2_clustering_header = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)

    }

    //
    // Create IGV session
    //
    IGV (
        ch_fasta,
        UCSC_BEDGRAPHTOBIGWIG.out.bigwig.collect{it[1]}.ifEmpty([]),
        MACS2_CALLPEAK.out.peak.collect{it[1]}.ifEmpty([]),
        MACS2_CONSENSUS.out.bed.collect{it[1]}.ifEmpty([]),
        params.modules['ucsc_bedgraphtobigwig'],
        params.modules['macs2_callpeak'],
        params.modules['macs2_consensus']
    )

    //
    // MODULE: Pipeline reporting
    //
    // ch_software_versions
    //     .map { it -> if (it) [ it.baseName, it ] }
    //     .groupTuple()
    //     .map { it[1][0] }
    //     .flatten()
    //     .collect()
    //     .set { ch_software_versions }

    // GET_SOFTWARE_VERSIONS (
    //     ch_software_versions.map { it }.collect()  //,
    //     // params.modules['get_software_versions'] // TODO check if its needed, I think that addParams does the trick
    // )
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile()
    )


    OUTPUT_DOCUMENTATION (
        ch_output_docs,
        ch_output_docs_images
    )

    //
    // MODULE: MultiQC
    //
    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowChipseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
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
            PHANTOMPEAKQUALTOOLS.out.spp.collect{it[1]}.ifEmpty([]),
            MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.nsc.collect{it[1]}.ifEmpty([]),
            MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.rsc.collect{it[1]}.ifEmpty([]),
            MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.correlation.collect{it[1]}.ifEmpty([]),

            MULTIQC_CUSTOM_PEAKS.out.count.collect{it[1]}.ifEmpty([]),
            MULTIQC_CUSTOM_PEAKS.out.frip.collect{it[1]}.ifEmpty([]),
            PLOT_HOMER_ANNOTATEPEAKS.out.tsv.collect().ifEmpty([]),
            SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]}.ifEmpty([]),
            // path ('macs/consensus/*') from ch_macs_consensus_deseq_mqc.collect().ifEmpty([])
        )
        multiqc_report       = MULTIQC.out.report.toList()
    }
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
