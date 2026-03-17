process PLOT_MACS3_QC {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b258138808975f2d51a13c5e4333cc9756eefbe0ce0e39ec90d5e0d69556bb39/data':
        'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-complexheatmap_r-base_r-ggplot2_pruned:1182f03e8fce2848' }"

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
