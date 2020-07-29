/*
 * Output Markdown documentation to HTML
 */
process OUTPUT_DOCUMENTATION {
    publishDir "${params.outdir}/${options.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (options.publish_results == "none") null
                      else filename }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path output_docs
    path images
    val options

    output:
    path "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}
