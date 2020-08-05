// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process MULTIQC_CUSTOM_PEAKS {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    tuple val(meta), path(peak), path(frip)
    path peak_count_header
    path frip_score_header
    val options

    output:
    tuple val(meta), path("*.peak_count_mqc.tsv"), emit: count
    tuple val(meta), path("*.FRiP_mqc.tsv"), emit: frip

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    """
    cat $peak | wc -l | awk -v OFS='\t' '{ print "${prefix}", \$1 }' | cat $peak_count_header - > ${prefix}.peak_count_mqc.tsv
    cat $frip_score_header $frip > ${prefix}.FRiP_mqc.tsv
    """
}
