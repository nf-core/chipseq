def SOFTWARE = 'homer'
def VERSION = '4.11'

process HOMER_ANNOTATEPEAKS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${options.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (options.publish_results == "none") null
                      else if (filename.endsWith('.version.txt')) null
                      else filename }

    container "quay.io/biocontainers/homer:4.11--pl526h9a982cc_2"
    //container "https://depot.galaxyproject.org/singularity/homer:4.11--pl526h9a982cc_2"

    conda (params.conda ? "bioconda::homer=4.11" : null)

    input:
    tuple val(meta), path(peak)
    path fasta
    path gtf
    val options

    output:
    tuple val(meta), path("*annotatePeaks.txt"), emit: txt
    path "*.version.txt", emit: version

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    annotatePeaks.pl \\
        $peak \\
        $fasta \\
        $options.args \\
        -gtf $gtf \\
        -cpu $task.cpus \\
        > ${prefix}.annotatePeaks.txt

    echo $VERSION > ${SOFTWARE}.version.txt
    """
}
