/*
 * Picard MergeSamFiles
 */
process PICARD_MERGESAMFILES {
    tag "$name"
    label 'process_medium'
    publishDir "${params.outdir}/bwa/mergedLibrary", mode: params.publish_dir_mode

    input:
    tuple val(name), val(single_end), path(bams)

    output:
    tuple val(name), val(single_end), path("*.bam")

    script:
    prefix = "${name}.mLb.mkD"
    bam_files = bams.sort()
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard MergeSamFiles] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.toGiga()
    }
    if (bam_files.size() > 1) {
        """
        picard -Xmx${avail_mem}g MergeSamFiles \\
            ${'INPUT='+bam_files.join(' INPUT=')} \\
            OUTPUT=${name}.sorted.bam \\
            SORT_ORDER=coordinate \\
            VALIDATION_STRINGENCY=LENIENT \\
            TMP_DIR=tmp
        """
    } else {
        """
        ln -s ${bam_files[0]} ${name}.sorted.bam
        """
    }
}
