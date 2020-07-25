/*
 * Convert GTF file to BED format
 */
process GTF2BED {
    tag "$gtf"
    label 'process_low'
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path gtf
    val opts

    output:
    path '*.bed'

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    """
    gtf2bed $gtf > ${gtf.baseName}.bed
    """
}
