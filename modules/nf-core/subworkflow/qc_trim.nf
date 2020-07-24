/*
 * Read QC and trimming
 */

include { FASTQC     } from '../software/fastqc'
include { TRIMGALORE } from '../software/trimgalore'

workflow QC_TRIM {
    take:
    ch_reads         // channel: [ val(meta), [ reads ] ]
    skip_fastqc      // boolean: true/false
    skip_trimming    // boolean: true/false
    fastqc_opts      //     map: options for FastQC module
    trim_galore_opts //     map: options for TrimGalore! module

    main:
    fastqc_html = Channel.empty()
    fastqc_zip = Channel.empty()
    fastqc_version = Channel.empty()
    if (!skip_fastqc) {
        FASTQC(ch_reads, fastqc_opts).html.set { fastqc_html }
        fastqc_zip = FASTQC.out.zip
        fastqc_version = FASTQC.out.version
    }

    ch_trim_reads = ch_reads
    trim_html = Channel.empty()
    trim_zip = Channel.empty()
    trim_log = Channel.empty()
    trim_galore_version = Channel.empty()
    if (!skip_trimming) {
        TRIMGALORE(ch_reads, trim_galore_opts).reads.set { ch_trim_reads }
        trim_html = TRIMGALORE.out.html
        trim_zip = TRIMGALORE.out.zip
        trim_log = TRIMGALORE.out.log
        trim_galore_version = TRIMGALORE.out.version
    }

    emit:
    fastqc_html           // channel: [ val(meta), [ html ] ]
    fastqc_zip            // channel: [ val(meta), [ zip ] ]
    fastqc_version        //    path: *.version.txt

    reads = ch_trim_reads // channel: [ val(meta), [ reads ] ]
    trim_html             // channel: [ val(meta), [ html ] ]
    trim_zip              // channel: [ val(meta), [ zip ] ]
    trim_log              // channel: [ val(meta), [ txt ] ]
    trim_galore_version   //    path: *.version.txt
}
