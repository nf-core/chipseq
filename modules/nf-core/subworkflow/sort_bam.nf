/*
 * Sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { SAMTOOLS_SORT  } from '../software/samtools_sort'
include { SAMTOOLS_INDEX } from '../software/samtools_index'
include { BAM_STATS      } from './bam_stats'

workflow SORT_BAM {
    take:
    ch_bam        // channel: [ val(meta), [ bam ] ]
    samtools_opts //     map: options for SAMTools modules

    main:
    SAMTOOLS_SORT(ch_bam, samtools_opts)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam, samtools_opts)
    BAM_STATS(SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]), samtools_opts)

    emit:
    bam = SAMTOOLS_SORT.out.bam                  // channel: [ val(meta), [ bam ] ]
    bai = SAMTOOLS_INDEX.out.bai                 // channel: [ val(meta), [ bai ] ]
    stats = BAM_STATS.out.stats                  // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS.out.flagstat            // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS.out.idxstats            // channel: [ val(meta), [ idxstats ] ]
    samtools_version = SAMTOOLS_SORT.out.version //    path: *.version.txt
}
