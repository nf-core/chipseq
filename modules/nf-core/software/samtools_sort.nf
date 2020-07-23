process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'
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
    tuple val(meta), path("*.bam"), emit: bam
    path "*.version.txt", emit: version

    script:
    prefix = opts.suffix ? "${meta.id}${opts.suffix}" : "${meta.id}"
    """
    samtools sort $opts.args -@ $task.cpus -o ${prefix}.bam -T $prefix $bam
    samtools --version | sed -n "s/.*\\(v.*\$\\)/\\1/p" > samtools.version.txt
    """
}
