// Import generic module functions
include { initOptions; saveFiles } from './functions'

/*
 * Parse software version numbers
 */
process GET_SOFTWARE_VERSIONS {
    publishDir "${params.outdir}/${options.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename, options, "pipeline_info") }

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
