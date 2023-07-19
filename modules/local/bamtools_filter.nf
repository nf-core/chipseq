process BAMTOOLS_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bamtools=2.5.2 bioconda::samtools=1.15.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0560a8046fc82aa4338588eca29ff18edab2c5aa:5687a7da26983502d0a8a9a6b05ed727c740ddc4-0' :
        'biocontainers/mulled-v2-0560a8046fc82aa4338588eca29ff18edab2c5aa:5687a7da26983502d0a8a9a6b05ed727c740ddc4-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed
    path bamtools_filter_se_config
    path bamtools_filter_pe_config

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        bamtools: \$(echo \$(bamtools --version 2>&1) | sed 's/^.*bamtools //; s/Part .*\$//')
    END_VERSIONS
    """
}
