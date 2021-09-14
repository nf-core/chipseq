// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
// custom_runName = params.name // TODO remove
// if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
//     custom_runName = workflow.runName
// }

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
    path ('trimgalore/*')
    path ('trimgalore/fastqc/*')

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
    path "*_data", emit: data
    // TODO add the line below:
    // path "*_plots"             , optional:true, emit: plots

    script:
    def software      = getSoftwareName(task.process)
    // def rtitle        = custom_runName ? "--title \"$custom_runName\"" : ''
    // def rfilename     = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    def custom_config = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f $options.args $custom_config .
    """
}
