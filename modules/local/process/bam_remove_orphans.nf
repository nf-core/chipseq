// Import generic module functions
include { initOptions; saveFiles } from './functions'

/*
 * Remove orphan reads from paired-end BAM file
 */
process BAM_REMOVE_ORPHANS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${options.publish_dir}${options.publish_by_id ? "/${meta.id}" : ''}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename, options, task.process.toLowerCase()) }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    tuple val(meta), path(bam)
    val options

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def ioptions = initOptions(options, task.process.toLowerCase())
    prefix = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    if (!meta.single_end) {
        """
        samtools sort -n -@ $task.cpus -o ${prefix}.name.sorted.bam -T ${prefix}.name.sorted $bam
        bampe_rm_orphan.py ${prefix}.name.sorted.bam ${prefix}.bam $ioptions.args
        """
    } else {
        """
        ln -s $bam ${prefix}.bam
        """
    }
}
