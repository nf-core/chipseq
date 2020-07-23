container "quay.io/biocontainers/macs2:2.2.7.1--py37h516909a_0"
//container "https://depot.galaxyproject.org/singularity/macs2:2.2.7.1--py37h516909a_0"
//echo \$(macs2 --version 2>&1) > macs2.version.txt


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
