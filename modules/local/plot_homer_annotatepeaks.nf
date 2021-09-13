// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Aggregated QC plots for peak-to-gene annotation
 */
process PLOT_HOMER_ANNOTATEPEAKS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path annos
    path mqc_header
    val suffix

    output:
    path '*.txt', emit: txt
    path '*.pdf', emit: pdf
    path '*.tsv', emit: tsv

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    """
    plot_homer_annotatepeaks.r \\
        -i ${annos.join(',')} \\
        -s ${annos.join(',').replaceAll("${suffix}","")} \\
        $options.args

    find ./ -type f -name "*.txt" -exec cat {} \\; | cat $mqc_header - > annotatepeaks.summary_mqc.tsv
    """
}
