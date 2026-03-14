process BAMTOOLS_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b935ba4da9d5bcaa5b036493218c4d76f4a018751fc6db47f144e88c16128eb5/data'
        : 'community.wave.seqera.io/library/bamtools_samtools:d0efa2e3de1e3441'}"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed
    path bamtools_filter_se_config
    path bamtools_filter_pe_config

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools
    tuple val("${task.process}"), val('bamtools'), eval("bamtools --version | sed '2!d;s/bamtools //g'"), topic: versions, emit: versions_bamtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def blacklist = bed ? "-L $bed" : ''
    def config = meta.single_end ? bamtools_filter_se_config : bamtools_filter_pe_config
    """
    samtools view \\
        $args \\
        $blacklist \\
        -b $bam \\
        | bamtools filter \\
            -out ${prefix}.bam \\
            -script $config
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
