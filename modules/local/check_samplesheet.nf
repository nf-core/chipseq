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

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def get_samplesheet_paths(LinkedHashMap row, boolean single_end) {
    def meta = [:]
    meta.id = row.sample_id
    meta.single_end = single_end

    def array = []
    if (single_end) {
        array = [ meta, [ file(row.fastq_1, checkIfExists: true) ] ]
    } else {
        array = [ meta, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ]
    }
    return array
}

// Function to get list of [sample, control, antibody, replicatesExist?, multipleGroups?]
def get_samplesheet_design(LinkedHashMap row, boolean single_end) {
    def meta = [:]
    meta.id = row.sample_id
    meta.single_end = single_end
    meta.antibody = row.antibody
    meta.control_id = row.control_id
    meta.reps_exist = row.replicatesExist.toBoolean()
    meta.groups_exist = row.multipleGroups.toBoolean()

    //return [ row.sample_id, row.control_id, row.antibody, row.replicatesExist.toBoolean(), row.multipleGroups.toBoolean() ]
    return meta
}
