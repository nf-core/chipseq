// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

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
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.conda ? "${baseDir}/environment.yml" : null)
    // conda (params.enable_conda ? "conda-forge::python=3.7.6 bioconda::samtools=1.9" : null) //TODO Create updated mulled container
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/mulled-v2-1a35167f7a491c7086c13835aaa74b39f1f43979:b8bbb1dc1cb75db128173603ca7a050729cdc526-0"
    // } else {
    //     container "quay.io/biocontainers/mulled-v2-1a35167f7a491c7086c13835aaa74b39f1f43979:b8bbb1dc1cb75db128173603ca7a050729cdc526-0"
    // } // TODO: Needs pysam library that is not in the mulled container

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (!meta.single_end) {
        """
        samtools sort -n -@ $task.cpus -o ${prefix}.name.sorted.bam -T ${prefix}.name.sorted $bam
        bampe_rm_orphan.py ${prefix}.name.sorted.bam ${prefix}.bam $options.args
        """
    } else {
        """
        ln -s $bam ${prefix}.bam
        """
    }
}
