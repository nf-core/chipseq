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

    container "quay.io/biocontainers/perl:5.26.2"
    //container "https://depot.galaxyproject.org/singularity/perl:5.26.2"

    conda (params.conda ? "conda-forge::perl=5.26.2" : null)

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
