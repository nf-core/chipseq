process HOMER_ANNOTATEPEAKS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else if (filename.endsWith('.version.txt')) null
                      else filename }

    container "quay.io/biocontainers/homer:4.11--pl526h9a982cc_2"
    //container "https://depot.galaxyproject.org/singularity/homer:4.11--pl526h9a982cc_2"

    conda (params.conda ? "bioconda::homer=4.11" : null)

    input:
    tuple val(meta), path(bed)
    path fasta
    path gtf
    val opts

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "*.version.txt", emit: version

    script:
    prefix = opts.suffix ? "${meta.id}${opts.suffix}" : "${meta.id}"

    """
    annotatePeaks.pl \\
        $opts.args \\
        $bed \\
        $fasta \\
        -gtf $gtf \\
        -cpu $task.cpus \\
        > ${prefix}.annotatePeaks.txt

    touch homer.version.txt
    """
}
