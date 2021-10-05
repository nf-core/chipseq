/*
 * Picard MarkDuplicates, sort, index BAM file and run samtools stats, flagstat and idxstats
 */
params.markduplicates_options = [:]
params.samtools_index_options = [:]
params.samtools_stats_options = [:]

include { PICARD_MARKDUPLICATES } from '../../modules/nf-core/modules/picard/markduplicates/main' addParams( options: params.markduplicates_options )
include { SAMTOOLS_INDEX        } from '../../modules/nf-core/modules/samtools/index/main'        addParams( options: params.samtools_index_options )
include { BAM_STATS_SAMTOOLS    } from './bam_stats_samtools'                                     addParams( options: params.samtools_stats_options )

workflow MARK_DUPLICATES_PICARD {
    take:
    ch_bam                 // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    PICARD_MARKDUPLICATES(ch_bam)
    SAMTOOLS_INDEX(PICARD_MARKDUPLICATES.out.bam)
    BAM_STATS_SAMTOOLS(PICARD_MARKDUPLICATES.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]))

    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first(),
                    SAMTOOLS_INDEX.out.versions.first(),
                    BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam               = PICARD_MARKDUPLICATES.out.bam      // channel: [ val(meta), [ bam ] ]
    metrics           = PICARD_MARKDUPLICATES.out.metrics  // channel: [ val(meta), [ metrics ] ]
    bai               = SAMTOOLS_INDEX.out.bai             // channel: [ val(meta), [ bai ] ]
    stats             = BAM_STATS_SAMTOOLS.out.stats       // channel: [ val(meta), [ stats ] ]
    flagstat          = BAM_STATS_SAMTOOLS.out.flagstat    // channel: [ val(meta), [ flagstat ] ]
    idxstats          = BAM_STATS_SAMTOOLS.out.idxstats    // channel: [ val(meta), [ idxstats ] ]

    versions          = ch_versions                        //    path: versions.yml
}
