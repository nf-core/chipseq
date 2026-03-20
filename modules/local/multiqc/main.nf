process MULTIQC {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34e733a9ae16a27e80fe00f863ea1479c96416017f24a907996126283e7ecd4d/data'
        : 'community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b'}"

    input:
    path workflow_summary
    path multiqc_config
    path mqc_custom_config
    path logo

    path "fastqc/*"
    path "trimgalore/fastqc/*"
    path "trimgalore/*"

    path "alignment/library/*"
    path "alignment/library/*"
    path "alignment/library/*"

    path "alignment/merged_library/unfiltered/*"
    path "alignment/merged_library/unfiltered/*"
    path "alignment/merged_library/unfiltered/*"
    path "alignment/merged_library/unfiltered/picard_metrics/*"

    path "alignment/merged_library/filtered/*"
    path "alignment/merged_library/filtered/*"
    path "alignment/merged_library/filtered/*"
    path "alignment/merged_library/filtered/picard_metrics/*"

    path "preseq/*"

    path "deeptools/*"
    path "deeptools/*"

    path "phantompeakqualtools/*"
    path "phantompeakqualtools/*"
    path "phantompeakqualtools/*"
    path "phantompeakqualtools/*"

    path "macs3/peaks/*"
    path "macs3/peaks/*"
    path "macs3/annotation/*"
    path "macs3/featurecounts/*"

    path "deseq2/*"
    path "deseq2/*"

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional: true, emit: plots
    path "versions.yml"        , emit: versions
    // From the nf-core/moduel: MultiQC should not push its versions to the `versions` topic. Its input depends on the versions topic to be resolved thus outputting to the topic will let the pipeline hang forever
    tuple val("${task.process}"), val('multiqc'), eval('multiqc --version | sed "s/.* //g"'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def custom_config = params.multiqc_config ? "--config ${mqc_custom_config}" : ''
    """
    multiqc \\
        -f \\
        ${args} \\
        ${custom_config} \\
        .
    """

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html
    touch multiqc_data/multiqc.log
    """
}
