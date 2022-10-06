/*
 * Map reads, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { BOWTIE2_ALIGN     } from '../../modules/nf-core/bowtie2/align/main'
include { BAM_SORT_SAMTOOLS } from './bam_sort_samtools'

workflow ALIGN_BOWTIE2 {
    take:
    reads          // channel: [ val(meta), [ reads ] ]
    index          //    path: /path/to/index
    save_unaligned // boolean: true/false

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with BWA
    //
    BOWTIE2_ALIGN(reads, index, save_unaligned, false)
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_SAMTOOLS(BOWTIE2_ALIGN.out.bam)
    ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions.first())

    emit:
    bam               = BAM_SORT_SAMTOOLS.out.bam               // channel: [ val(meta), [ bam ] ]
    bai               = BAM_SORT_SAMTOOLS.out.bai               // channel: [ val(meta), [ bai ] ]
    stats             = BAM_SORT_SAMTOOLS.out.stats             // channel: [ val(meta), [ stats ] ]
    flagstat          = BAM_SORT_SAMTOOLS.out.flagstat          // channel: [ val(meta), [ flagstat ] ]
    idxstats          = BAM_SORT_SAMTOOLS.out.idxstats          // channel: [ val(meta), [ idxstats ] ]

    versions          = ch_versions                             //    path: versions.yml
}
