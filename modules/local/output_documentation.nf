// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Output Markdown documentation to HTML
 */
process OUTPUT_DOCUMENTATION {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }

    conda (params.conda ? "${baseDir}/environment.yml" : null) // TODO update with pointers to singularity and docker container

    input:
    path output_docs
    path images

    output:
    path "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}
