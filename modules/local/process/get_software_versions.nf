/*
 * Parse software version numbers
 */
process GET_SOFTWARE_VERSIONS {
  publishDir "${params.outdir}/${options.publish_dir}",
      mode: params.publish_dir_mode,
      saveAs: { filename ->
                    if (options.publish_results == "none") null
                    else if (filename.endsWith('.yaml')) null
                    else filename }

    input:
    path versions
    val options

    output:
    path "software_versions.csv", emit: csv
    path 'software_versions_mqc.yaml', emit: yaml

    script:
    """
    echo $workflow.manifest.version > pipeline.version.txt
    echo $workflow.nextflow.version > nextflow.version.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}
