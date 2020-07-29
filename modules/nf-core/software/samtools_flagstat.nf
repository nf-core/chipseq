def SOFTWARE = 'samtools'

process SAMTOOLS_FLAGSTAT {
    tag "$meta.id"
    publishDir "${params.outdir}/${options.publish_dir}${options.publish_by_id ? "/${meta.id}" : ''}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (options.publish_results == "none") null
                      else if (filename.endsWith('.version.txt')) null
                      else filename }

    container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    //container " https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"

    conda (params.conda ? "bioconda::samtools=1.10" : null)

    input:
    tuple val(meta), path(bam), path(bai)
    val options

    output:
    tuple val(meta), path("*.flagstat"), emit: flagstat
    path "*.version.txt", emit: version

    script:
    """
    samtools flagstat $bam > ${bam}.flagstat
    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${SOFTWARE}.version.txt
    """
}
