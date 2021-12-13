process BEDTOOLS_GENOMECOV {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0':
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(bam), path(flagstat)

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    tuple val(meta), path("*.txt")     , emit: scale_factor
    path "versions.yml"                , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def pe     = meta.single_end ? '' : '-pc'
    def extend = (meta.single_end && params.fragment_size > 0) ? "-fs ${params.fragment_size}" : ''
    """
    SCALE_FACTOR=\$(grep '[0-9] mapped (' $flagstat | awk '{print 1000000/\$1}')
    echo \$SCALE_FACTOR > ${prefix}.scale_factor.txt

    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -scale \$SCALE_FACTOR \\
        $pe \\
        $extend \\
        | sort -T '.' -k1,1 -k2,2n > ${prefix}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
