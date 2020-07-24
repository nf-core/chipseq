/*
 * Check input samplesheet and get read channels
 */

include { CHECK_SAMPLESHEET;
          get_samplesheet_paths;
          get_samplesheet_design } from '../process/check_samplesheet'

workflow CHECK_INPUT {
    take:
    ch_input               //   file: /path/to/samplesheet.csv
    seq_center             // string: sequencing center for read group
    check_samplesheet_opts //    map: options for TrimGalore! module

    main:
    CHECK_SAMPLESHEET (ch_input, check_samplesheet_opts)
        .reads
        .splitCsv(header:true, sep:',')
        .map { get_samplesheet_paths(it, params.single_end, seq_center) }
        .set { ch_reads }

    CHECK_SAMPLESHEET
        .out
        .controls
        .splitCsv(header:true, sep:',')
        .map { get_samplesheet_design(it, params.single_end) }
        .set { ch_design }

    emit:
    reads = ch_reads   // channel: [ val(meta), [ reads ] ]
    design = ch_design // channel: [ val(meta), [ reads ] ]
}
