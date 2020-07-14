/*
 * Index BAM file
 */
process SAMTOOLS_INDEX {
    tag "$name"
    label 'process_medium'
    publishDir path: "${params.outdir}", mode: params.publish_dir_mode
    
    input:
    tuple val(name), val(single_end), path(bam)

    output:
    tuple val(name), val(single_end), path('*.bai')

    script:
    prefix = "${name}.Lb"
    """
    samtools index ${prefix}.sorted.bam
    """
}
