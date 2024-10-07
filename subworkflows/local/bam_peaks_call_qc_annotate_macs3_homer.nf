//
// Call peaks with MACS3, annotate with HOMER and perform downstream QC
//

include { MACS3_CALLPEAK           } from '../../modules/nf-core/macs3/callpeak/main'
include { HOMER_ANNOTATEPEAKS      } from '../../modules/nf-core/homer/annotatepeaks/main'

include { FRIP_SCORE               } from '../../modules/local/frip_score'
include { MULTIQC_CUSTOM_PEAKS     } from '../../modules/local/multiqc_custom_peaks'
include { PLOT_MACS3_QC            } from '../../modules/local/plot_macs3_qc'
include { PLOT_HOMER_ANNOTATEPEAKS } from '../../modules/local/plot_homer_annotatepeaks'

workflow BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER {
    take:
    ch_bam                            // channel: [ val(meta), [ ip_bam ], [ control_bam ] ]
    ch_fasta                          // channel: [ fasta ]
    ch_gtf                            // channel: [ gtf ]
    macs_gsize                        // integer: value for --macs_gsize parameter
    annotate_peaks_suffix             //  string: suffix for input HOMER annotate peaks files to be trimmed off
    ch_peak_count_header_multiqc      // channel: [ header_file ]
    ch_frip_score_multiqc             // channel: [ header_file ]
    ch_peak_annotation_header_multiqc // channel: [ header_file ]
    is_narrow_peak                    // boolean: true/false
    skip_peak_annotation              // boolean: true/false
    skip_peak_qc                      // boolean: true/false

    main:

    ch_versions = Channel.empty()

    //
    // Call peaks with MACS3
    //
    MACS3_CALLPEAK (
        ch_bam,
        macs_gsize
    )
    ch_versions = ch_versions.mix(MACS3_CALLPEAK.out.versions.first())

    //
    // Filter out samples with 0 MACS3 peaks called
    //
    MACS3_CALLPEAK
        .out
        .peak
        .filter {
            meta, peaks ->
                peaks.size() > 0
        }
        .set { ch_macs3_peaks }

    // Create channels: [ meta, ip_bam, peaks ]
    ch_bam
        .join(ch_macs3_peaks, by: [0])
        .map {
            meta, ip_bam, control_bam, peaks ->
                [ meta, ip_bam, peaks ]
        }
        .set { ch_bam_peaks }

    //
    // Calculate FRiP score
    //
    FRIP_SCORE (
        ch_bam_peaks
    )
    ch_versions = ch_versions.mix(FRIP_SCORE.out.versions.first())

    // Create channels: [ meta, peaks, frip ]
    ch_bam_peaks
        .join(FRIP_SCORE.out.txt, by: [0])
        .map {
            meta, ip_bam, peaks, frip ->
                [ meta, peaks, frip ]
        }
        .set { ch_bam_peak_frip }

    //
    // FRiP score custom content for MultiQC
    //
    MULTIQC_CUSTOM_PEAKS (
        ch_bam_peak_frip,
        ch_peak_count_header_multiqc,
        ch_frip_score_multiqc
    )
    ch_versions = ch_versions.mix(MULTIQC_CUSTOM_PEAKS.out.versions.first())

    ch_homer_annotatepeaks          = Channel.empty()
    ch_plot_macs3_qc_txt            = Channel.empty()
    ch_plot_macs3_qc_pdf            = Channel.empty()
    ch_plot_homer_annotatepeaks_txt = Channel.empty()
    ch_plot_homer_annotatepeaks_pdf = Channel.empty()
    ch_plot_homer_annotatepeaks_tsv = Channel.empty()
    if (!skip_peak_annotation) {
        //
        // Annotate peaks with HOMER
        //
        HOMER_ANNOTATEPEAKS (
            ch_macs3_peaks,
            ch_fasta,
            ch_gtf
        )
        ch_homer_annotatepeaks = HOMER_ANNOTATEPEAKS.out.txt
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS.out.versions.first())

        if (!skip_peak_qc) {
            //
            // MACS3 QC plots with R
            //
            PLOT_MACS3_QC (
                ch_macs3_peaks.collect{it[1]},
                is_narrow_peak
            )
            ch_plot_macs3_qc_txt = PLOT_MACS3_QC.out.txt
            ch_plot_macs3_qc_pdf = PLOT_MACS3_QC.out.pdf
            ch_versions = ch_versions.mix(PLOT_MACS3_QC.out.versions)

            //
            // Peak annotation QC plots with R
            //
            PLOT_HOMER_ANNOTATEPEAKS (
                HOMER_ANNOTATEPEAKS.out.txt.collect{it[1]},
                ch_peak_annotation_header_multiqc,
                annotate_peaks_suffix
            )
            ch_plot_homer_annotatepeaks_txt = PLOT_HOMER_ANNOTATEPEAKS.out.txt
            ch_plot_homer_annotatepeaks_pdf = PLOT_HOMER_ANNOTATEPEAKS.out.pdf
            ch_plot_homer_annotatepeaks_tsv = PLOT_HOMER_ANNOTATEPEAKS.out.tsv
            ch_versions = ch_versions.mix(PLOT_HOMER_ANNOTATEPEAKS.out.versions)
        }
    }

    emit:
    peaks                        = ch_macs3_peaks                   // channel: [ val(meta), [ peaks ] ]
    xls                          = MACS3_CALLPEAK.out.xls           // channel: [ val(meta), [ xls ] ]
    gapped_peaks                 = MACS3_CALLPEAK.out.gapped        // channel: [ val(meta), [ gapped_peak ] ]
    bed                          = MACS3_CALLPEAK.out.bed           // channel: [ val(meta), [ bed ] ]
    bedgraph                     = MACS3_CALLPEAK.out.bdg           // channel: [ val(meta), [ bedgraph ] ]

    frip_txt                     = FRIP_SCORE.out.txt               // channel: [ val(meta), [ txt ] ]

    frip_multiqc                 = MULTIQC_CUSTOM_PEAKS.out.frip    // channel: [ val(meta), [ frip ] ]
    peak_count_multiqc           = MULTIQC_CUSTOM_PEAKS.out.count   // channel: [ val(meta), [ counts ] ]

    homer_annotatepeaks          = ch_homer_annotatepeaks           // channel: [ val(meta), [ txt ] ]

    plot_macs3_qc_txt            = ch_plot_macs3_qc_txt             // channel: [ txt ]
    plot_macs3_qc_pdf            = ch_plot_macs3_qc_pdf             // channel: [ pdf ]

    plot_homer_annotatepeaks_txt = ch_plot_homer_annotatepeaks_txt  // channel: [ txt ]
    plot_homer_annotatepeaks_pdf = ch_plot_homer_annotatepeaks_pdf  // channel: [ pdf ]
    plot_homer_annotatepeaks_tsv = ch_plot_homer_annotatepeaks_tsv  // channel: [ tsv ]

    versions                     = ch_versions                      // channel: [ versions.yml ]
}
