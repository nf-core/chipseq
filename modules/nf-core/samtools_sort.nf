/*
 * Sort BAM file
 */
process SAMTOOLS_SORT {
    tag "$name"
    label 'process_medium'
    publishDir path: "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(name), val(single_end), path(bam)
    
    output:
    tuple val(name), val(single_end), path('*.bam')

    script:
    prefix = "${name}.Lb"
    """
    samtools sort -@ $task.cpus -o ${prefix}.sorted.bam -T $name $bam
    """
}
