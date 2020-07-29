def SOFTWARE = 'deeptools'

process DEEPTOOLS_PLOTHEATMAP {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/${options.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (options.publish_results == "none") null
                      else if (filename.endsWith('.version.txt')) null
                      else filename }

    container "quay.io/biocontainers/deeptools:3.4.3--py_0"
    //container "https://depot.galaxyproject.org/singularity/deeptools:3.4.3--py_0"

    conda (params.conda ? "bioconda::deeptools=3.4.3" : null)

    input:
    tuple val(meta), path(matrix)
    val options

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.tab"), emit: table
    path "*.version.txt", emit: version

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    plotHeatmap \\
        $options.args \\
        --matrixFile $matrix \\
        --outFileName ${prefix}.plotHeatmap.pdf \\
        --outFileNameMatrix ${prefix}.plotHeatmap.mat.tab

    plotHeatmap --version | sed -e "s/plotHeatmap //g" > ${SOFTWARE}.version.txt
    """
}
