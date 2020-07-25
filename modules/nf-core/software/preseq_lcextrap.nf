process PRESEQ_LCEXTRAP {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else if (filename.endsWith('.version.txt')) null
                      else filename }

    container "quay.io/biocontainers/preseq:2.0.3--hf53bd2b_3"
    //container "https://depot.galaxyproject.org/singularity/preseq:2.0.3--hf53bd2b_3"

    conda (params.conda ? "bioconda::preseq=2.0.3" : null)

    input:
    tuple val(meta), path(bam)
    val opts

    output:
    tuple val(meta), path("*.ccurve.txt"), emit: ccurve
    tuple val(meta), path("*.log"), emit: log
    path "*.version.txt", emit: version

    script:
    prefix = opts.suffix ? "${meta.id}${opts.suffix}" : "${meta.id}"
    pe = meta.single_end ? '' : '-pe'
    """
    preseq \\
        lc_extrap \\
        $opts.args \\
        $pe \\
        -output ${prefix}.ccurve.txt \\
        $bam
    cp .command.err ${prefix}.command.log
    preseq &> preseq.version.txt
    """
}
