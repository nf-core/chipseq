/*
 * Reformat design file, check validitiy and create IP vs control mappings
 */
process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    path samplesheet

    output:
    path 'samplesheet_reads.csv', emit: reads
    path 'samplesheet_controls.csv', emit: controls

    script:  // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    """
    check_samplesheet.py $samplesheet samplesheet_reads.csv samplesheet_controls.csv
    """
}

// Function to get list of [ sample, single_end?, [ fastq_1, fastq_2 ] ]
def get_samplesheet_paths(LinkedHashMap row, boolean single_end) {
    def sample = row.sample_id
    def fastq_1 = row.fastq_1
    def fastq_2 = row.fastq_2

    def array = []
    if (single_end) {
        array = [ sample, single_end, [ file(fastq_1, checkIfExists: true) ] ]
    } else {
        array = [ sample, single_end, [ file(fastq_1, checkIfExists: true), file(fastq_2, checkIfExists: true) ] ]
    }
    return array
}

// Function to get list of [sample, control, antibody, replicatesExist?, multipleGroups?]
def get_samplesheet_design(LinkedHashMap row) {
    return [ row.sample_id, row.control_id, row.antibody, row.replicatesExist.toBoolean(), row.multipleGroups.toBoolean() ]
}
