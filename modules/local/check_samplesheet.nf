/*
 * Reformat design file, check validitiy and create IP vs control mappings
 */
process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    path samplesheet

    output:
    path 'design_reads.csv'
    path 'design_controls.csv'

    script:  // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    """
    check_samplesheet.py $design design_reads.csv design_controls.csv
    """
}

// // Function to get list of [ sample, single_end?, [ fastq_1, fastq_2 ] ]
// def check_samplesheet_paths(LinkedHashMap samplesheet) {
//     def sample = samplesheet.sample
//     def single_end = samplesheet.single_end.toBoolean()
//     def fastq_1 = samplesheet.fastq_1
//     def fastq_2 = samplesheet.fastq_2
//
//     def array = []
//     if (single_end) {
//         array = [ sample, single_end, [ file(fastq_1, checkIfExists: true) ] ]
//     } else {
//         array = [ sample, single_end, [ file(fastq_1, checkIfExists: true), file(fastq_2, checkIfExists: true) ] ]
//     }
//     return array
// }
