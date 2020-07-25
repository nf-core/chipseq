conda (params.conda ? "${baseDir}/environment.yml" : null)
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
