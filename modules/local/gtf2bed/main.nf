process GTF2BED {
    tag "$gtf"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2':
        'biocontainers/perl:5.26.2' }"

    input:
    path gtf

    output:
    path '*.bed'       , emit: bed
    tuple val("${task.process}"), val('perl'), eval("perl -V:version | sed \"s/version='//; s/';//\""), topic:versions, emit: versions_perl

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    """
    gtf2bed \\
        $gtf \\
        > ${gtf.baseName}.bed
    """
}
