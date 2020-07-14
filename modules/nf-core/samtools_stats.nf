/*
 * Index BAM file
 */
process SAMTOOLS_STATS {
    tag "$name"
    label 'process_medium'
    publishDir path: "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(name), val(single_end), path(bam)

    output:
    tuple val(name), val(single_end), path('*.stats'), emit: stats
    tuple val(name), val(single_end), path('*.flagstat'), emit: flagstat
    tuple val(name), val(single_end), path('*.idxstats'), emit: idxstats

    script:
    prefix = "${name}.Lb"
    """
    samtools stats $bam > ${bam}.stats
    samtools flagstat $bam > ${bam}.flagstat
    samtools idxstats $bam > ${bam}.idxstats
    """
}
