// Import generic module functions
include { initOptions; saveFiles } from './functions'

def SOFTWARE = 'preseq'

process PRESEQ_LCEXTRAP {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}/${options.publish_dir}${options.publish_by_id ? "/${meta.id}" : ''}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename, options, SOFTWARE) }

    container "quay.io/biocontainers/preseq:2.0.3--hf53bd2b_3"
    //container "https://depot.galaxyproject.org/singularity/preseq:2.0.3--hf53bd2b_3"

    conda (params.conda ? "bioconda::preseq=2.0.3" : null)

    input:
    tuple val(meta), path(bam)
    val options

    output:
    tuple val(meta), path("*.ccurve.txt"), emit: ccurve
    tuple val(meta), path("*.log"), emit: log
    path "*.version.txt", emit: version

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    pe = meta.single_end ? '' : '-pe'
    """
    preseq \\
        lc_extrap \\
        $options.args \\
        $pe \\
        -output ${prefix}.ccurve.txt \\
        $bam
    cp .command.err ${prefix}.command.log

    echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//' > ${SOFTWARE}.version.txt
    """
}
