/*
 * Filter BAM file
 */
process FILTER_BAM {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else if (filename.endsWith('.version.txt')) null
                    else filename }

    input:
    tuple val(meta), path(bam), path(bai)
    path bed
    path config
    val opts

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.version.txt", emit: version

    script:
    prefix = opts.suffix ? "${meta.id}${opts.suffix}" : "${meta.id}"
    filter_params = meta.single_end ? '-F 0x004' : '-F 0x004 -F 0x0008 -f 0x001'
    dup_params = params.keep_dups ? '' : '-F 0x0400'
    multimap_params = params.keep_multi_map ? '' : '-q 1'
    blacklist_params = params.blacklist ? "-L $bed" : ''
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
    echo \$(bamtools --version 2>&1) > bamtools.version.txt
    """
}
