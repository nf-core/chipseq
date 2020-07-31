process MULTIQC_CUSTOM_PEAKS {
    tag "$meta.id"
    publishDir "${params.outdir}/${options.publish_dir}${options.publish_by_id ? "/${meta.id}" : ''}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (options.publish_results == "none") null
                      else if (filename.endsWith('.version.txt')) null
                      else filename }

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
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    cat $peak | wc -l | awk -v OFS='\t' '{ print "${prefix}", \$1 }' | cat $peak_count_header - > ${prefix}.peak_count_mqc.tsv
    cat $frip_score_header $frip > ${prefix}.FRiP_mqc.tsv
    """
}
