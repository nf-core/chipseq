// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MULTIQC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::multiqc=1.10.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/multiqc:1.10.1--pyhdfd78af_1"
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--pyhdfd78af_1"
    }

    input:
    path multiqc_config
    path mqc_custom_config
    path software_versions
    path workflow_summary

    path ('fastqc/*')
    path ('trimgalore/fastqc/*')
    path ('trimgalore/*')

    path ('alignment/library/*')
    path ('alignment/library/*')
    path ('alignment/library/*')

    path ('alignment/mergedLibrary/unfiltered/*')
    path ('alignment/mergedLibrary/unfiltered/*')
    path ('alignment/mergedLibrary/unfiltered/*')
    path ('alignment/mergedLibrary/unfiltered/picard_metrics/*')

    path ('alignment/mergedLibrary/filtered/*')
    path ('alignment/mergedLibrary/filtered/*')
    path ('alignment/mergedLibrary/filtered/*')
    path ('alignment/mergedLibrary/filtered/picard_metrics/*')

    path ('preseq/*')
    path ('deeptools/*')
    path ('deeptools/*')
    path ('phantompeakqualtools/*')
    path ('phantompeakqualtools/*')
    path ('phantompeakqualtools/*')
    path ('phantompeakqualtools/*')

    path ('macs2/peaks/*')
    path ('macs2/peaks/*')
    path ('macs2/annotation/*')

    path ('featurecounts/*')
    // path ('macs/consensus/*') from ch_macs_consensus_deseq_mqc.collect().ifEmpty([])

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    script:
    def custom_config = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc \\
        -f \\
        $options.args \\
        $custom_config \\
        .

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
