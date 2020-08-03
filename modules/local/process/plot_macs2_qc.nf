// Import generic module functions
include { initOptions; saveFiles } from './functions'

/*
 * Aggregated QC plots for peaks
 */
process PLOT_MACS2_QC {
    label 'process_medium'
    publishDir "${params.outdir}/${options.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename, options, task.process.tokenize('_')[1].toLowerCase()) }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path peaks
    val options

    output:
    path '*.txt', emit: txt
    path '*.pdf', emit: pdf

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def software = task.process.tokenize('_')[1].toLowerCase()
    def ioptions = initOptions(options, software)
    peak_type = params.narrow_peak ? 'narrowPeak' : 'broadPeak'
    """
    plot_macs2_qc.r \\
        -i ${peaks.join(',')} \\
        -s ${peaks.join(',').replaceAll("_peaks.${peak_type}","")} \\
        $ioptions.args
    """
}
