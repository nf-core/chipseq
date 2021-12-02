process MULTIQC_CUSTOM_PEAKS {

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img':
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(peak), path(frip)
    path peak_count_header
    path frip_score_header

    output:
    tuple val(meta), path("*.peak_count_mqc.tsv"), emit: count
    tuple val(meta), path("*.FRiP_mqc.tsv"), emit: frip

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $peak | wc -l | awk -v OFS='\t' '{ print "${prefix}", \$1 }' | cat $peak_count_header - > ${prefix}.peak_count_mqc.tsv
    cat $frip_score_header $frip > ${prefix}.FRiP_mqc.tsv
    """
}
