def SOFTWARE = 'deeptools'

process DEEPTOOLS_COMPUTEMATRIX {
    tag "$meta.id"
    label 'process_high'
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
    tuple val(meta), path(bigwig)
    path bed
    val options

    output:
    tuple val(meta), path("*.mat.gz"), emit: matrix
    tuple val(meta), path("*.mat.tab"), emit: table
    path "*.version.txt", emit: version

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    computeMatrix \\
        $options.args \\
        --regionsFileName $bed \\
        --scoreFileName $bigwig \\
        --outFileName ${prefix}.computeMatrix.mat.gz \\
        --outFileNameMatrix ${prefix}.computeMatrix.vals.mat.tab \\
        --numberOfProcessors $task.cpus

    computeMatrix --version | sed -e "s/computeMatrix //g" > ${SOFTWARE}.version.txt
    """
}
