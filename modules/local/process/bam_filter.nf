// Import generic module functions
include { initOptions; saveFiles } from './functions'

/*
 * Filter BAM file
 */
process BAM_FILTER {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename=filename, options=options, publish_dir=task.process.toLowerCase(), publish_id=meta.id) }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    tuple val(meta), path(bam), path(bai)
    path bed
    path bamtools_filter_se_config
    path bamtools_filter_pe_config
    val options

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.version.txt", emit: version

    script:
    def ioptions = initOptions(options)
    prefix = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    filter_params = meta.single_end ? '-F 0x004' : '-F 0x004 -F 0x0008 -f 0x001'
    dup_params = params.keep_dups ? '' : '-F 0x0400'
    multimap_params = params.keep_multi_map ? '' : '-q 1'
    blacklist_params = params.blacklist ? "-L $bed" : ''
    config = meta.single_end ? bamtools_filter_se_config : bamtools_filter_pe_config
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
    echo \$(bamtools --version 2>&1) | sed 's/^.*bamtools //; s/Part .*\$//' > bamtools.version.txt
    """
}
