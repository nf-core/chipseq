process ANNOTATE_BOOLEAN_PEAKS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'docker.io/library/ubuntu:20.04' }"

    input:
    tuple val(meta), path(boolean_txt), path(homer_peaks)

    output:
    path '*.boolean.annotatePeaks.txt', emit: annotate_peaks_txt
    tuple val("${task.process}"), val('sed'), eval("sed --version 2>&1 | sed '1!d;s/^.*) //'"), topic: versions, emit: versions_sed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cut -f2- ${homer_peaks} | awk 'NR==1; NR > 1 {print \$0 | "sort -T '.' -k1,1 -k2,2n"}' | cut -f6- > tmp.txt
    paste $boolean_txt tmp.txt > ${prefix}.boolean.annotatePeaks.txt
    """
}
