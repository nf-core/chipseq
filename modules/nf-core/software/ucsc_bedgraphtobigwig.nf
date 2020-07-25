def SOFTWARE = 'ucsc'
def VERSION = '377'

process UCSC_BEDRAPHTOBIGWIG {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else if (filename.endsWith('.version.txt')) null
                      else filename }

    container "quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1"
    //container "https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:377--h446ed27_1"

    conda (params.conda ? "bioconda::ucsc-bedgraphtobigwig=377" : null)

    input:
    tuple val(meta), path(bedgraph)
    path sizes
    val opts

    output:
    tuple val(meta), path("*.bigWig"), emit: bigwig

    script:
    prefix = opts.suffix ? "${meta.id}${opts.suffix}" : "${meta.id}"
    """
    bedGraphToBigWig $bedgraph $sizes ${prefix}.bigWig
    echo $VERSION > ${SOFTWARE}.version.txt
    """
}
