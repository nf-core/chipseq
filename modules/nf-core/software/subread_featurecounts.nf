def SOFTWARE = 'subread'

process SUBREAD_FEATURECOUNTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${options.publish_dir}${options.publish_by_id ? "/${meta.id}" : ''}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (options.publish_results == "none") null
                      else if (filename.endsWith('.version.txt')) null
                      else filename }

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
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    pe = meta.single_end ? '' : '-p'
    """
    featureCounts \\
        $options.args \\
        $pe \\
        -T $task.cpus \\
        -a $annotation \\
        -o ${prefix}.featureCounts.txt \\
        ${bams.join(' ')}

    echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g" > ${SOFTWARE}.version.txt
    """
}
