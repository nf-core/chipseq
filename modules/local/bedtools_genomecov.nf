process BEDTOOLS_GENOMECOV {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else if (filename.endsWith('.version.txt')) null
                    else filename }

    container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    //container "https://depot.galaxyproject.org/singularity/bedtools:2.29.2--hc088bd4_0"

    conda (params.conda ? "${moduleDir}/environment.yml" : null)

    input:
    tuple val(meta), path(bam), path(flagstat)
    val opts

    output:
    tuple val(meta), path("*.bigWig"), emit: bigwig
    tuple val(meta), path("*.txt"), emit: scale_factor
    path "*.version.txt", emit: version

    script:
    prefix = opts.suffix ? "${meta.id}${opts.suffix}" : "${meta.id}"
    pe = meta.single_end ? '' : '-pc'
    extend = (params.single_end && params.fragment_size > 0) ? "-fs ${params.fragment_size}" : ''
    """
    SCALE_FACTOR=\$(grep 'mapped (' $flagstat | awk '{print 1000000/\$1}')
    echo \$SCALE_FACTOR > ${prefix}.scale_factor.txt

    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -scale \$SCALE_FACTOR \\
        $pe \\
        $extend \\
        | sort -T '.' -k1,1 -k2,2n >  ${prefix}.bedGraph

    bedtools --version > bedtools.version.txt
    """
}
//     tuple val(name), path(bam), path(flagstat) from ch_rm_orphan_bam_bigwig.join(ch_rm_orphan_flagstat_bigwig, by: [0])
//     path sizes from ch_genome_sizes_bigwig.collect()
// bedGraphToBigWig ${prefix}.bedGraph $sizes ${prefix}.bigWig
//find * -type f -name "*.bigWig" -exec echo -e "bwa/mergedLibrary/bigwig/"{}"\\t0,0,178" \\; > ${name}.bigWig.igv.txt
