// Import generic module functions
include { initOptions; saveFiles } from './functions'

def SOFTWARE = 'samtools'

process SAMTOOLS_IDXSTATS {
    tag "$meta.id"
    publishDir "${params.outdir}/${options.publish_dir}${options.publish_by_id ? "/${meta.id}" : ''}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename, options, SOFTWARE) }

    container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    //container " https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"

    conda (params.conda ? "bioconda::samtools=1.10" : null)

    input:
    tuple val(meta), path(bam), path(bai)
    val options

    output:
    tuple val(meta), path("*.idxstats"), emit: idxstats
    path "*.version.txt", emit: version

    script:
    """
    samtools idxstats $bam > ${bam}.idxstats
    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${SOFTWARE}.version.txt
    """
}
