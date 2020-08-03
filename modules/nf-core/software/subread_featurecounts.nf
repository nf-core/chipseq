// Import generic module functions
include { initOptions; saveFiles } from './functions'

process SUBREAD_FEATURECOUNTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${options.publish_dir}${options.publish_by_id ? "/${meta.id}" : ''}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename, options, task.process.tokenize('_')[0].toLowerCase()) }
        
    container "quay.io/biocontainers/subread:2.0.1--hed695b0_0"
    //container "https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0"

    conda (params.conda ? "bioconda::subread=2.0.1" : null)

    input:
    tuple val(meta), path(bams), path(annotation)
    val options

    output:
    tuple val(meta), path("*featureCounts.txt"), emit: txt
    tuple val(meta), path("*featureCounts.txt.summary"), emit: summary
    path "*.version.txt", emit: version

    script:
    def software = task.process.tokenize('_')[0].toLowerCase()
    def ioptions = initOptions(options, software)
    prefix = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    pe = meta.single_end ? '' : '-p'
    """
    featureCounts \\
        $ioptions.args \\
        $pe \\
        -T $task.cpus \\
        -a $annotation \\
        -o ${prefix}.featureCounts.txt \\
        ${bams.join(' ')}

    echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g" > ${software}.version.txt
    """
}
