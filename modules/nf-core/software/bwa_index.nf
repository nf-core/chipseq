def SOFTWARE = 'bwa'

// Function to initialise default values and to generate a Groovy Map of module options
def init_options (Map args, String publish_dir) {
    def Map options = [:]
    options.args            = args.args ?: ''
    options.args2           = args.args2 ?: ''
    options.suffix          = args.suffix ?: ''
    options.publish_dir     = args.publish_dir ? args.publish_dir : publish_dir
    options.publish_files   = args.publish_files ?: null
    options.publish_by_id   = args.publish_by_id ?: false
    return options
}

// Function to publish module results
  // if publish_files == null           - all files are published
  // if publish_files == Map [:]        - no files are published
  // if publish_files == Map [ext:path] - Only files that end with "ext" are published to "path" appended to output directory
def publish_options ( filename, publish_files ) {
    if (publish_files instanceof Map) {
        if (!publish_files.isEmpty()) {
            if (!filename.endsWith('.version.txt')) {
                filename
                // Loop through extensions and paths here
            }
        }
    }
    else {
        if (!filename.endsWith('.version.txt')) filename
    }
}

process BWA_INDEX {
    tag "$fasta"
    label 'process_high'
    //if (module_options.publish_files) {
        publishDir "${params.outdir}/${options.publish_dir}",
            mode: params.publish_dir_mode,
            saveAs: { filename -> publish_options( filename, module_options.publish_files ) }
    //}

    container "biocontainers/bwa:v0.7.17_cv1"
    //container "https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7"

    conda (params.conda ? "bioconda::bwa=0.7.17" : null)

    //echo true
    //cache false

    input:
    path fasta
    val options

    output:
    path "${fasta}.*", emit: index
    path "*.version.txt", emit: version

    script:
    def module_options = init_options(options, SOFTWARE)
    """
    echo $module_options
    bwa index $options.args $fasta
    echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//' > ${SOFTWARE}.version.txt
    """
}
