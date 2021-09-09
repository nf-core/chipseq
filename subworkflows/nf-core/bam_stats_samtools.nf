/*
 * Run SAMtools stats, flagstat and idxstats
 */

include { SAMTOOLS_STATS    } from '../../modules/nf-core/modules/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../modules/nf-core/modules/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/modules/samtools/flagstat/main'

workflow BAM_STATS_SAMTOOLS {
    take:
    ch_bam_bai       // channel: [ val(meta), [ bam ], [bai] ]
    samtools_options //     map: options for SAMTools modules

    main:
    SAMTOOLS_STATS(ch_bam_bai, samtools_options)
    SAMTOOLS_FLAGSTAT(ch_bam_bai, samtools_options)
    SAMTOOLS_IDXSTATS(ch_bam_bai, samtools_options)

    emit:
    stats = SAMTOOLS_STATS.out.stats              // channel: [ val(meta), [ stats ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat     // channel: [ val(meta), [ flagstat ] ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats     // channel: [ val(meta), [ idxstats ] ]
    samtools_version = SAMTOOLS_STATS.out.version //    path: *.version.txt
}
