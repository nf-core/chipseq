// Import generic module functions
include { initOptions; saveFiles } from './functions'

process DEEPTOOLS_COMPUTEMATRIX {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/${options.publish_dir}${options.publish_by_id ? "/${meta.id}" : ''}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename, options, task.process.tokenize('_')[0].toLowerCase()) }

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
    def software = task.process.tokenize('_')[0].toLowerCase()
    def ioptions = initOptions(options, software)
    prefix = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    """
    computeMatrix \\
        $ioptions.args \\
        --regionsFileName $bed \\
        --scoreFileName $bigwig \\
        --outFileName ${prefix}.computeMatrix.mat.gz \\
        --outFileNameMatrix ${prefix}.computeMatrix.vals.mat.tab \\
        --numberOfProcessors $task.cpus

    computeMatrix --version | sed -e "s/computeMatrix //g" > ${software}.version.txt
    """
}
