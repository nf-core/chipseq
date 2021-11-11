// Import generic module functions
include { initOptions; saveFiles; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Remove orphan reads from paired-end BAM file
 */
process BAM_REMOVE_ORPHANS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:task.process.toLowerCase(), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::pysam=0.16.0.1 bioconda::samtools=1.13" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-57736af1eb98c01010848572c9fec9fff6ffaafd:980645707aa512b63b32c44c35888ac1dab126dd-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-57736af1eb98c01010848572c9fec9fff6ffaafd:980645707aa512b63b32c44c35888ac1dab126dd-0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    path "versions.yml"                   , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (!meta.single_end) {
        """
        samtools sort -n -@ $task.cpus -o ${prefix}.name.sorted.bam -T ${prefix}.name.sorted $bam
        bampe_rm_orphan.py ${prefix}.name.sorted.bam ${prefix}.bam $options.args

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        ln -s $bam ${prefix}.bam

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}
