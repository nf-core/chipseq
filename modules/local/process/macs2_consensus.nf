// conda (params.conda ? "${baseDir}/environment.yml" : null)
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
