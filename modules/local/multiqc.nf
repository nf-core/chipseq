process MULTIQC {
    label 'process_medium'
    conda "bioconda::multiqc=1.25.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0'
        : 'biocontainers/multiqc:1.25.1--pyhdfd78af_0'}"

    input:
    path workflow_summary
    path multiqc_config
    path mqc_custom_config
    path logo

    path "fastqc/*"
    path "trimgalore/fastqc/*"
    path "trimgalore/*"

    path "alignment/library/samtools_stats/*"
    path "alignment/library/samtools_flagstat/*"
    path "alignment/library/samtools_idxstats/*"

    path "alignment/merged_library/unfiltered/samtools_stats/*"
    path "alignment/merged_library/unfiltered/samtools_flagstat/*"
    path "alignment/merged_library/unfiltered/samtools_idxstats/*"
    path "alignment/merged_library/unfiltered/picard_metrics/*"

    path "alignment/merged_library/filtered/samtools_stats/*"
    path "alignment/merged_library/filtered/samtools_flagstat/*"
    path "alignment/merged_library/filtered/samtools_idxstats/*"
    path "alignment/merged_library/filtered/picard_metrics/*"

    path "preseq/*"

    path "deeptools/plotprofile/*"
    path "deeptools/plotfingerprint/*"

    path "phantompeakqualtools/spp/*"
    path "phantompeakqualtools/nsc/*"
    path "phantompeakqualtools/rsc/*"
    path "phantompeakqualtools/correlation/*"

    path "macs3/peaks/frip/*"
    path "macs3/peaks/count/*"
    path "macs3/annotation/*"
    path "macs3/featurecounts/*"

    path "deseq2/pca/*"
    path "deseq2/clustering/*"

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional: true, emit: plots
    path "versions.yml"        , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html
    touch multiqc_data/multiqc.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
