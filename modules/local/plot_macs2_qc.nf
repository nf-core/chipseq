// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Aggregated QC plots for peaks
 */
process PLOT_MACS2_QC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path peaks

    output:
    path '*.txt', emit: txt
    path '*.pdf', emit: pdf

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def peak_type = params.narrow_peak ? 'narrowPeak' : 'broadPeak'
    """
    plot_macs2_qc.r \\
        -i ${peaks.join(',')} \\
        -s ${peaks.join(',').replaceAll("_peaks.${peak_type}","")} \\
        $options.args
    """
}
