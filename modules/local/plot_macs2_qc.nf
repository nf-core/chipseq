// Import generic module functions
include { initOptions; saveFiles; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Aggregated QC plots for peaks
 */
process PLOT_MACS2_QC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:task.process.toLowerCase(), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "${baseDir}/environment.yml" : null) //TODO Create updated mulled container

    input:
    path peaks

    output:
    path '*.txt'       , emit: txt
    path '*.pdf'       , emit: pdf
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def peak_type = params.narrow_peak ? 'narrowPeak' : 'broadPeak'
    """
    plot_macs2_qc.r \\
        -i ${peaks.join(',')} \\
        -s ${peaks.join(',').replaceAll("_peaks.${peak_type}","")} \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        python: \$(python --version | sed 's/Python //g')
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
