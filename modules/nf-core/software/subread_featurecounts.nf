container "quay.io/biocontainers/subread:2.0.1--hed695b0_0"
//container "https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0"
//echo \$(featureCounts -v 2>&1) > featurecounts.version.txt
conda (params.conda ? "bioconda::subread=2.0.1" : null)

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
