process SAMTOOLS_FLAGSTAT {
    tag "$meta.id"
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else if (filename.endsWith('.version.txt')) null
                      else filename }

    container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    //container " https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"

    conda (params.conda ? "bioconda::samtools=1.10" : null)

    input:
    tuple val(meta), path(bam), path(bai)
    val opts

    output:
    tuple val(meta), path("*.flagstat"), emit: flagstat
    path "*.version.txt", emit: version

    script:
    """
    samtools flagstat $bam > ${bam}.flagstat
    samtools --version | sed -n "s/.*\\(v.*\$\\)/\\1/p" > samtools.version.txt
    """
}
