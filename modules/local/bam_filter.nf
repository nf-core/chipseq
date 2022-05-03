/*
 * Filter BAM file
 */
process BAM_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bamtools=2.5.2 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0560a8046fc82aa4338588eca29ff18edab2c5aa:5687a7da26983502d0a8a9a6b05ed727c740ddc4-0' :
        'quay.io/biocontainers/mulled-v2-0560a8046fc82aa4338588eca29ff18edab2c5aa:5687a7da26983502d0a8a9a6b05ed727c740ddc4-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed
    path bamtools_filter_se_config
    path bamtools_filter_pe_config

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def filter_params    = meta.single_end ? '-F 0x004' : '-F 0x004 -F 0x0008 -f 0x001'
    def dup_params       = params.keep_dups ? '' : '-F 0x0400'
    def multimap_params  = params.keep_multi_map ? '' : '-q 1'
    def blacklist_params = params.blacklist ? "-L $bed" : ''
    def config           = meta.single_end ? bamtools_filter_se_config : bamtools_filter_pe_config
    """
    samtools view \\
        $filter_params \\
        $dup_params \\
        $multimap_params \\
        $blacklist_params \\
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
