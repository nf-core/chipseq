/*
 * Remove orphan reads from paired-end BAM file
 */
process REMOVE_BAM_ORPHANS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else filename }

    conda (params.conda ? "${baseDir}/environment.yml" : null)
    
    input:
    tuple val(meta), path(bam)
    val opts

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    prefix = opts.suffix ? "${meta.id}${opts.suffix}" : "${meta.id}"
    if (!meta.single_end) {
        """
        samtools sort -n -@ $task.cpus -o ${prefix}.name.sorted.bam -T ${prefix}.name.sorted $bam
        bampe_rm_orphan.py ${prefix}.name.sorted.bam ${prefix}.bam $opts.args
        """
    } else {
        """
        ln -s $bam ${prefix}.bam
        """
    }
}
