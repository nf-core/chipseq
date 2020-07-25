/*
 * Picard MarkDuplicates, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { PICARD_MARKDUPLICATES } from '../software/picard_markduplicates'
include { SAMTOOLS_INDEX        } from '../software/samtools_index'
include { BAM_STATS             } from './bam_stats'

workflow MARK_DUPLICATES {
    take:
    ch_bam              // channel: [ val(meta), [ bam ] ]
    markduplicates_opts //     map: options for picard MarkDuplicates module
    samtools_opts       //     map: options for SAMTools modules

    main:
    PICARD_MARKDUPLICATES(ch_bam, markduplicates_opts)
    SAMTOOLS_INDEX(PICARD_MARKDUPLICATES.out.bam, samtools_opts)
    BAM_STATS(PICARD_MARKDUPLICATES.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]), samtools_opts)

    emit:
    bam = PICARD_MARKDUPLICATES.out.bam                // channel: [ val(meta), [ bam ] ]
    metrics = PICARD_MARKDUPLICATES.out.metrics        // channel: [ val(meta), [ metrics ] ]
    bai = SAMTOOLS_INDEX.out.bai                       // channel: [ val(meta), [ bai ] ]
    stats = BAM_STATS.out.stats                        // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS.out.flagstat                  // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS.out.idxstats                  // channel: [ val(meta), [ idxstats ] ]
    picard_version = PICARD_MARKDUPLICATES.out.version //    path: *.version.txt
    samtools_version = SAMTOOLS_INDEX.out.version      //    path: *.version.txt
}
