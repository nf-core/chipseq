// Import generic module functions
include { initOptions; saveFiles } from './functions'

/*
 * Aggregated QC plots for peak-to-gene annotation
 */
process PLOT_HOMER_ANNOTATEPEAKS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename=filename, options=options, publish_dir=task.process.toLowerCase(), publish_id='') }
        
    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path annos
    path mqc_header
    val suffix
    val options

    output:
    path '*.txt', emit: txt
    path '*.pdf', emit: pdf
    path '*.tsv', emit: tsv

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def ioptions = initOptions(options)
    """
    plot_homer_annotatepeaks.r \\
        -i ${annos.join(',')} \\
        -s ${annos.join(',').replaceAll("${suffix}","")} \\
        $ioptions.args

    find ./ -type f -name "*.txt" -exec cat {} \\; | cat $mqc_header - > annotatepeaks.summary_mqc.tsv
    """
}
