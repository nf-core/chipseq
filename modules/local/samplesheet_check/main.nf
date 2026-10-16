process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/62/622d8944750bc95bb56b4c3ed5c2b827e677c14073d48a5231e0f2bec0718add/data' :
        'community.wave.seqera.io/library/python:3.12.12--74abbf3898230efd' }"

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args = task.ext.args ?: ''
    """
    check_samplesheet.py \\
        $samplesheet \\
        $args
    """
}
