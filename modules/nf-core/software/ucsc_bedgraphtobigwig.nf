def SOFTWARE = 'ucsc'
def VERSION = '377'

process UCSC_BEDRAPHTOBIGWIG {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${options.publish_dir}${options.publish_by_id ? "/${meta.id}" : ''}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (options.publish_results == "none") null
                      else if (filename.endsWith('.version.txt')) null
                      else filename }

    container "quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1"
    //container "https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:377--h446ed27_1"

    conda (params.conda ? "bioconda::ucsc-bedgraphtobigwig=377" : null)

    input:
    tuple val(meta), path(bedgraph)
    path sizes
    val options

    output:
    tuple val(meta), path("*.bigWig"), emit: bigwig
    path "*.version.txt", emit: version

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    bedGraphToBigWig $bedgraph $sizes ${prefix}.bigWig
    echo $VERSION > ${SOFTWARE}.version.txt
    """
}
