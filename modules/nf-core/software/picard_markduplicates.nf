process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else if (filename.endsWith('.version.txt')) null
                    else filename }

    container "quay.io/biocontainers/picard:2.23.2--0"
    //container "https://depot.galaxyproject.org/singularity/picard:2.23.2--0"

    conda (params.conda ? "bioconda::picard=2.23.2" : null)
    
    input:
    tuple val(meta), path(bam)
    val opts

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.metrics.txt", emit: metrics
    path "*.version.txt", emit: version

    script:
    prefix = opts.suffix ? "${meta.id}${opts.suffix}" : "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.toGiga()
    }
    """
    picard \\
        -Xmx${avail_mem}g \\
        MarkDuplicates \\
        $opts.args \\
        INPUT=$bam \\
        OUTPUT=${prefix}.bam \\
        METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt
    picard MarkDuplicates --version &> picard.version.txt || true
    """
}
