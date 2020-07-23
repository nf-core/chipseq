process SAMTOOLS_INDEX {
    tag "$meta.id"
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else if (filename.endsWith('.version.txt')) null
                    else filename }

    container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    //container " https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"

    conda (params.conda ? "${moduleDir}/environment.yml" : null)

    input:
    tuple val(meta), path(bam)
    val opts

    output:
    tuple val(meta), path("*.bai"), emit: bai
    path "*.version.txt", emit: version

    script:
    """
    samtools index $bam
    samtools --version | sed -n "s/.*\\(v.*\$\\)/\\1/p" > samtools.version.txt
    """
}
