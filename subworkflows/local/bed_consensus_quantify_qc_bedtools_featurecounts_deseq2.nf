//
// Call consensus peaks with BEDTools and custom scripts, annotate with HOMER, quantify with featureCounts and QC with DESeq2
//

include { HOMER_ANNOTATEPEAKS    } from '../../modules/nf-core/homer/annotatepeaks/main'
include { SUBREAD_FEATURECOUNTS  } from '../../modules/nf-core/subread/featurecounts/main'

include { MACS3_CONSENSUS        } from '../../modules/local/macs3_consensus'
include { ANNOTATE_BOOLEAN_PEAKS } from '../../modules/local/annotate_boolean_peaks'
include { DESEQ2_QC              } from '../../modules/local/deseq2_qc'

workflow BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2 {
    take:
    ch_peaks                            // channel: [ val(meta), [ peaks ] ]
    ch_bams                             // channel: [ val(meta), [ ip_bams ] ]
    ch_fasta                            // channel: [ fasta ]
    ch_gtf                              // channel: [ gtf ]
    ch_deseq2_pca_header_multiqc        // channel: [ header_file ]
    ch_deseq2_clustering_header_multiqc // channel: [ header_file ]
    is_narrow_peak                      // boolean: true/false
    skip_peak_annotation                // boolean: true/false
    skip_deseq2_qc                      // boolean: true/false

    main:

    ch_versions = Channel.empty()

    // Create channels: [ meta , [ peaks ] ]
    // Where meta = [ id:antibody, multiple_groups:true/false, replicates_exist:true/false ]
    ch_peaks
        .map {
            meta, peak ->
                [ meta.antibody, meta.id - ~/_REP\d+$/, peak ]
        }
        .groupTuple()
        .map {
            antibody, groups, peaks ->
                [
                    antibody,
                    groups.groupBy().collectEntries { [(it.key) : it.value.size()] },
                    peaks
                ]
        }
        .map {
            antibody, groups, peaks ->
                def meta_new = [:]
                meta_new.id = antibody
                meta_new.multiple_groups = groups.size() > 1
                meta_new.replicates_exist = groups.max { it.value }.value > 1
                [ meta_new, peaks ]
        }
        .set { ch_antibody_peaks }

    //
    // Generate consensus peaks across samples
    //
    MACS3_CONSENSUS (
        ch_antibody_peaks,
        is_narrow_peak
    )
    ch_versions = ch_versions.mix(MACS3_CONSENSUS.out.versions)

    //
    // Annotate consensus peaks
    //
    if (!skip_peak_annotation) {
        HOMER_ANNOTATEPEAKS (
            MACS3_CONSENSUS.out.bed,
            ch_fasta,
            ch_gtf
        )
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS.out.versions)

        //
        // MODULE: Add boolean fields to annotated consensus peaks to aid filtering
        //
        ANNOTATE_BOOLEAN_PEAKS (
            MACS3_CONSENSUS.out.boolean_txt.join(HOMER_ANNOTATEPEAKS.out.txt, by: [0]),
        )
        ch_versions = ch_versions.mix(ANNOTATE_BOOLEAN_PEAKS.out.versions)
    }

    // Create channels: [ meta, [ ip_bams ], saf ]
    MACS3_CONSENSUS
        .out
        .saf
        .map {
            meta, saf ->
                [ meta.id, meta, saf ]
        }
        .join(ch_bams)
        .map {
            antibody, meta, saf, bams, meta_single_end ->
                [ meta + meta_single_end, bams.flatten().sort(), saf ]
        }
        .set { ch_bam_saf }

    //
    // Quantify peaks across samples with featureCounts
    //
    SUBREAD_FEATURECOUNTS (
        ch_bam_saf
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions)

    //
    // Generate QC plots with DESeq2
    //
    ch_deseq2_qc_pdf           = Channel.empty()
    ch_deseq2_qc_rdata         = Channel.empty()
    ch_deseq2_qc_rds           = Channel.empty()
    ch_deseq2_qc_pca_txt       = Channel.empty()
    ch_deseq2_qc_pca_multiqc   = Channel.empty()
    ch_deseq2_qc_dists_txt     = Channel.empty()
    ch_deseq2_qc_dists_multiqc = Channel.empty()
    ch_deseq2_qc_log           = Channel.empty()
    ch_deseq2_qc_size_factors  = Channel.empty()
    if (!skip_deseq2_qc) {
        DESEQ2_QC (
            SUBREAD_FEATURECOUNTS.out.counts,
            ch_deseq2_pca_header_multiqc,
            ch_deseq2_clustering_header_multiqc
        )
        ch_deseq2_qc_pdf           = DESEQ2_QC.out.pdf
        ch_deseq2_qc_rdata         = DESEQ2_QC.out.rdata
        ch_deseq2_qc_rds           = DESEQ2_QC.out.rds
        ch_deseq2_qc_pca_txt       = DESEQ2_QC.out.pca_txt
        ch_deseq2_qc_pca_multiqc   = DESEQ2_QC.out.pca_multiqc
        ch_deseq2_qc_dists_txt     = DESEQ2_QC.out.dists_txt
        ch_deseq2_qc_dists_multiqc = DESEQ2_QC.out.dists_multiqc
        ch_deseq2_qc_log           = DESEQ2_QC.out.log
        ch_deseq2_qc_size_factors  = DESEQ2_QC.out.size_factors
        ch_versions = ch_versions.mix(DESEQ2_QC.out.versions)
    }

    emit:
    consensus_bed           = MACS3_CONSENSUS.out.bed           // channel: [ bed ]
    consensus_saf           = MACS3_CONSENSUS.out.saf           // channel: [ saf ]
    consensus_pdf           = MACS3_CONSENSUS.out.pdf           // channel: [ pdf ]
    consensus_txt           = MACS3_CONSENSUS.out.txt           // channel: [ pdf ]
    consensus_boolean_txt   = MACS3_CONSENSUS.out.boolean_txt   // channel: [ txt ]
    consensus_intersect_txt = MACS3_CONSENSUS.out.intersect_txt // channel: [ txt ]

    featurecounts_txt       = SUBREAD_FEATURECOUNTS.out.counts  // channel: [ txt ]
    featurecounts_summary   = SUBREAD_FEATURECOUNTS.out.summary // channel: [ txt ]

    deseq2_qc_pdf           = ch_deseq2_qc_pdf                  // channel: [ pdf ]
    deseq2_qc_rdata         = ch_deseq2_qc_rdata                // channel: [ rdata ]
    deseq2_qc_rds           = ch_deseq2_qc_rds                  // channel: [ rds ]
    deseq2_qc_pca_txt       = ch_deseq2_qc_pca_txt              // channel: [ txt ]
    deseq2_qc_pca_multiqc   = ch_deseq2_qc_pca_multiqc          // channel: [ txt ]
    deseq2_qc_dists_txt     = ch_deseq2_qc_dists_txt            // channel: [ txt ]
    deseq2_qc_dists_multiqc = ch_deseq2_qc_dists_multiqc        // channel: [ txt ]
    deseq2_qc_log           = ch_deseq2_qc_log                  // channel: [ txt ]
    deseq2_qc_size_factors  = ch_deseq2_qc_size_factors         // channel: [ txt ]

    versions                = ch_versions                       // channel: [ versions.yml ]
}
