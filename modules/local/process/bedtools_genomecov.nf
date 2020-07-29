def SOFTWARE = 'bedtools'

process BEDTOOLS_GENOMECOV {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${options.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (options.publish_results == "none") null
                      else if (filename.endsWith('.version.txt')) null
                      else filename }

    container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    //container "https://depot.galaxyproject.org/singularity/bedtools:2.29.2--hc088bd4_0"

    conda (params.conda ? "bioconda::bedtools=2.29.2" : null)

    input:
    tuple val(meta), path(bam), path(flagstat)
    val options

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    tuple val(meta), path("*.txt"), emit: scale_factor
    path "*.version.txt", emit: version

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    pe = meta.single_end ? '' : '-pc'
    extend = (meta.single_end && params.fragment_size > 0) ? "-fs ${params.fragment_size}" : ''
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
        | sort -T '.' -k1,1 -k2,2n > ${prefix}.bedGraph

    bedtools --version | sed -e "s/bedtools v//g" > ${SOFTWARE}.version.txt
    """
}
