process PICARD_MERGESAMFILES {
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
    tuple val(meta), path(bams)
    val opts

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.version.txt", emit: version

    script:
    prefix = opts.suffix ? "${meta.id}${opts.suffix}" : "${meta.id}"
    bam_files = bams.sort()
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard MergeSamFiles] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.toGiga()
    }
    if (bam_files.size() > 1) {
        """
        picard \\
            -Xmx${avail_mem}g \\
            MergeSamFiles \\
            $opts.args \\
            ${'INPUT='+bam_files.join(' INPUT=')} \\
            OUTPUT=${prefix}.bam
        picard MergeSamFiles --version &> picard.version.txt || true
        """
    } else {
        """
        ln -s ${bam_files[0]} ${prefix}.bam
        picard MergeSamFiles --version &> picard.version.txt || true
        """
    }
}
