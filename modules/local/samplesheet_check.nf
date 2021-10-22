//TODO Update with rnaseq 3.3 version
// Import generic module functions
include { saveFiles; getProcessName } from './functions'

/*
 * Reformat design file, check validitiy and create IP vs control mappings
 */
process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:"pipeline_info", publish_id:'') }

    conda (params.enable_conda ? "${baseDir}/environment.yml" : null)

    input:
    path samplesheet
    val options

    output:
    path '*.csv'

    script:  // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    """
    check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
