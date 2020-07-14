/*
 * Picard MarkDuplicates
 */
process PICARD_MARKDUPLICATES {
    tag "$name"
    label 'process_medium'
    publishDir "${params.outdir}/bwa/mergedLibrary", mode: params.publish_dir_mode

    input:
    tuple val(name), val(single_end), path(bam)

    output:
    tuple val(name), val(single_end), path("*.bam"), emit: bam
    tuple val(name), val(single_end), path("*.txt"), emit: metrics

    script:
    prefix = "${name}.mLb.mkD"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.toGiga()
    }
    """
    picard -Xmx${avail_mem}g MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${prefix}.sorted.bam \\
        ASSUME_SORTED=true \\
        REMOVE_DUPLICATES=false \\
        METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt \\
        VALIDATION_STRINGENCY=LENIENT \\
        TMP_DIR=tmp
    """
}
