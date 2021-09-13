/*
 * Read QC and trimming
 */


params.fastqc_options     = [:]
params.trimgalore_options = [:]

include { FASTQC     } from '../../modules/nf-core/modules/fastqc/main'     addParams( options: params.fastqc_options )
include { TRIMGALORE } from '../../modules/nf-core/modules/trimgalore/main' addParams( options: params.trimgalore_options )

workflow FASTQC_TRIMGALORE {
    take:
    ch_reads           // channel: [ val(meta), [ reads ] ]
    skip_fastqc        // boolean: true/false
    skip_trimming      // boolean: true/false

    main:
    fastqc_html = Channel.empty()
    fastqc_zip = Channel.empty()
    fastqc_version = Channel.empty()
    if (!skip_fastqc) {
        FASTQC(ch_reads).html.set { fastqc_html }
        fastqc_zip = FASTQC.out.zip
        fastqc_version = FASTQC.out.version
    }

    ch_trim_reads = ch_reads
    trim_html = Channel.empty()
    trim_zip = Channel.empty()
    trim_log = Channel.empty()
    trimgalore_version = Channel.empty()
    if (!skip_trimming) {
        TRIMGALORE(ch_reads).reads.set { ch_trim_reads }
        trim_html = TRIMGALORE.out.html
        trim_zip = TRIMGALORE.out.zip
        trim_log = TRIMGALORE.out.log
        trimgalore_version = TRIMGALORE.out.version
    }

    emit:
    fastqc_html           // channel: [ val(meta), [ html ] ]
    fastqc_zip            // channel: [ val(meta), [ zip ] ]
    fastqc_version        //    path: *.version.txt

    reads = ch_trim_reads // channel: [ val(meta), [ reads ] ]
    trim_html             // channel: [ val(meta), [ html ] ]
    trim_zip              // channel: [ val(meta), [ zip ] ]
    trim_log              // channel: [ val(meta), [ txt ] ]
    trimgalore_version    //    path: *.version.txt
}
