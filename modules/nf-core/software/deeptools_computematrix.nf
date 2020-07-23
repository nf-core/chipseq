process DEEPTOOLS_COMPUTEMATRIX {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else if (filename.endsWith('.version.txt')) null
                    else filename }

    container "quay.io/biocontainers/deeptools:3.4.3--py_0"
    //container "https://depot.galaxyproject.org/singularity/deeptools:3.4.3--py_0"

    conda (params.conda ? "${moduleDir}/environment.yml" : null)

    input:
    tuple val(meta), path(bigwig)
    path bed
    val opts

    output:
    tuple val(meta), path("*.mat.gz"), emit: matrix
    tuple val(meta), path("*.mat.tab"), emit: table
    path "*.version.txt", emit: version

    script:
    prefix = opts.suffix ? "${meta.id}${opts.suffix}" : "${meta.id}"
    """
    computeMatrix \\
        $opts.args \\
        --regionsFileName $bed \\
        --scoreFileName $bigwig \\
        --outFileName ${prefix}.computeMatrix.mat.gz \\
        --outFileNameMatrix ${prefix}.computeMatrix.vals.mat.tab \\
        --numberOfProcessors $task.cpus

    echo \$(computeMatrix --version 2>&1) > deeptools.version.txt || true
    """
}
