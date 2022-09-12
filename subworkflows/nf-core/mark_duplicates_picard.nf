/*
 * Picard MarkDuplicates, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { PICARD_MARKDUPLICATES } from '../../modules/nf-core/modules/picard/markduplicates/main'
include { SAMTOOLS_INDEX        } from '../../modules/nf-core/modules/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from './bam_stats_samtools'

workflow MARK_DUPLICATES_PICARD {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    //
    // Picard MarkDuplicates
    //
    PICARD_MARKDUPLICATES(bam)
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX(PICARD_MARKDUPLICATES.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    BAM_STATS_SAMTOOLS(PICARD_MARKDUPLICATES.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]))
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = PICARD_MARKDUPLICATES.out.bam     // channel: [ val(meta), [ bam ] ]
    metrics  = PICARD_MARKDUPLICATES.out.metrics // channel: [ val(meta), [ metrics ] ]

    bai      = SAMTOOLS_INDEX.out.bai            // channel: [ val(meta), [ bai ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                       // channel: [ versions.yml ]
}
