/*
 * Get chromosome sizes from a fasta file
 */
process GET_CHROM_SIZES {
    tag "$fasta"
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    conda (params.conda ? "${baseDir}/environment.yml" : null)
    
    input:
    path fasta
    val opts

    output:
    path '*.sizes', emit: sizes
    path '*.fai', emit: fai

    script:
    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
    """
}
