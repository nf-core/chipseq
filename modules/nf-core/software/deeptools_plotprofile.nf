process DEEPTOOLS_PLOTPROFILE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else if (filename.endsWith('.version.txt')) null
                    else filename }

    container "quay.io/biocontainers/deeptools:3.4.3--py_0"
    //container "https://depot.galaxyproject.org/singularity/deeptools:3.4.3--py_0"

    conda (params.conda ? "bioconda::deeptools=3.4.3" : null)
    
    input:
    tuple val(meta), path(matrix)
    val opts

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.tab"), emit: table
    path "*.version.txt", emit: version

    script:
    prefix = opts.suffix ? "${meta.id}${opts.suffix}" : "${meta.id}"
    """
    plotProfile \\
        $opts.args \\
        --matrixFile $matrix \\
        --outFileName ${prefix}.plotProfile.pdf \\
        --outFileNameData ${prefix}.plotProfile.tab

    echo \$(plotProfile --version 2>&1) > deeptools.version.txt || true
    """
}
