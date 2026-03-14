/*
 * Remove orphan reads from paired-end BAM file
 */
process BAM_REMOVE_ORPHANS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container 
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/87b06acf94e50ddd0d4ce705952ecb496a81b6665f7971cd7492f270d921cb6c/data'
        : 'community.wave.seqera.io/library/pysam_samtools:b9e3a5f6b6caee59'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    if (!meta.single_end) {
        """
        samtools sort -n -@ $task.cpus -o ${prefix}.name.sorted.bam -T ${prefix}.name.sorted $bam
        bampe_rm_orphan.py ${prefix}.name.sorted.bam ${prefix}.bam $args
        """
    } else {
        """
        ln -s $bam ${prefix}.bam
        """
    }

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
