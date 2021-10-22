//
// Check input samplesheet and get read channels
//

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

/*
 * Check input samplesheet and get read channels
 */
workflow INPUT_CHECK {
    take:
    ch_input                  //   file: /path/to/samplesheet.csv
    seq_center                // string: sequencing center for read group
    samplesheet_check_options //    map: options for check_samplesheet module

    main:
    SAMPLESHEET_CHECK (ch_input, samplesheet_check_options)
        .splitCsv(header:true, sep:',')
        .map { get_samplesheet_paths(it, seq_center) }
        .set { ch_reads }

    emit:
    reads = ch_reads // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def get_samplesheet_paths(LinkedHashMap row, String seq_center) {
    def meta = [:]
    meta.id = row.sample
    meta.single_end = row.single_end.toBoolean()
    meta.antibody = row.antibody
    meta.control = row.control

    def rg = "\'@RG\\tID:${meta.id}\\tSM:${meta.id.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:1\'"
    if (seq_center) {
        rg = "\'@RG\\tID:${meta.id}\\tSM:${meta.id.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:1\\tCN:${seq_center}\'"
    }
    meta.read_group = rg

    def array = []
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1, checkIfExists: true) ] ]
    } else {
        array = [ meta, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ]
    }
    return array
}
