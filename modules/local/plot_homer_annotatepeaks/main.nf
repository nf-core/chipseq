process PLOT_HOMER_ANNOTATEPEAKS {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/23/23e639aceaf0ba484c7d306cbd39c038117c12d35beccfa6e811c96d22bea799/data':
        'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-complexheatmap_r-base_r-ggplot2_pruned:d46617a59d25f733' }"

    input:
    path annos
    path mqc_header
    val suffix

    output:
    path '*.txt'       , emit: txt
    path '*.pdf'       , emit: pdf
    path '*.tsv'       , emit: tsv
    tuple val("${task.process}"), val('R'), eval('R --version | sed "1!d; s/.*version //; s/ .*//"'), topic: versions, emit: versions_r

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "annotatepeaks"
    """
    plot_homer_annotatepeaks.r \\
        -i ${annos.join(',')} \\
        -s ${annos.join(',').replaceAll("${suffix}","")} \\
        -p $prefix \\
        $args

    find ./ -type f -name "*summary.txt" -exec cat {} \\; | cat $mqc_header - > ${prefix}.summary_mqc.tsv
    """
}
