process PLOT_MACS3_QC {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/23/23e639aceaf0ba484c7d306cbd39c038117c12d35beccfa6e811c96d22bea799/data':
        'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-complexheatmap_r-base_r-ggplot2_pruned:d46617a59d25f733' }"

    input:
    path peaks
    val is_narrow_peak

    output:
    path '*.txt'       , emit: txt
    path '*.pdf'       , emit: pdf
    tuple val("${task.process}"), val('R'), eval('R --version | sed "1!d; s/.*version //; s/ .*//"'), topic: versions, emit: versions_r

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args      = task.ext.args ?: ''
    def peak_type = is_narrow_peak ? 'narrowPeak' : 'broadPeak'
    """
    plot_macs3_qc.r \\
        -i ${peaks.join(',')} \\
        -s ${peaks.join(',').replaceAll("_peaks.${peak_type}","")} \\
        $args
    """
}
