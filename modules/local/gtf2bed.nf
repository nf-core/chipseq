/*
 * Convert GTF to BED file
 */
process GTF2BED {
    tag "$gtf"
    //label 'process_low'
    publishDir "${params.outdir}/genome", mode: params.publish_dir_mode

    input:
    path gtf

    output:
    path '*.bed'

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    """
    gtf2bed $gtf > ${gtf.baseName}.bed
    """
}
