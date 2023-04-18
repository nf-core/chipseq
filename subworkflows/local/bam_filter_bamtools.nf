include { SAMTOOLS_SORT           } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX          } from '../../modules/nf-core/samtools/index/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'
include { BAM_STATS_SAMTOOLS      } from '../nf-core/bam_stats_samtools/main'

include { BAMTOOLS_FILTER         } from '../../modules/local/bamtools_filter'
include { BAM_REMOVE_ORPHANS      } from '../../modules/local/bam_remove_orphans'

workflow BAM_FILTER_BAMTOOLS {
    take:
    ch_bam_bai                   // channel: [ val(meta), [ bam ], [bai] ]
    ch_bed                       // channel: [ bed ]
    ch_fasta                     // channel: [ fasta ]
    ch_bamtools_filter_se_config // channel: [ config_file ]
    ch_bamtools_filter_pe_config // channel: [ config_file ]

    main:

    ch_versions = Channel.empty()

    //
    // Filter BAM file with BAMTools
    //
    BAMTOOLS_FILTER (
        ch_bam_bai,
        ch_bed,
        ch_bamtools_filter_se_config,
        ch_bamtools_filter_pe_config
    )
    ch_versions = ch_versions.mix(BAMTOOLS_FILTER.out.versions.first())

    BAMTOOLS_FILTER
        .out
        .bam
        .branch {
            meta, bam ->
                single_end: meta.single_end
                    return [ meta, bam ]
                paired_end: !meta.single_end
                    return [ meta, bam ]
        }
        .set { ch_bam }

    //
    // Index SE BAM file
    //
    SAMTOOLS_INDEX {
        ch_bam.single_end
    }
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // Run samtools stats, flagstat and idxstats on SE BAM
    //
    BAM_STATS_SAMTOOLS (
        ch_bam.single_end.join(SAMTOOLS_INDEX.out.bai),
        ch_fasta
    )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions.first())

    //
    // Name sort PE BAM before filtering with pysam
    //
    SAMTOOLS_SORT (
        ch_bam.paired_end
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    //
    // Remove orphan reads from PE BAM file
    //
    BAM_REMOVE_ORPHANS (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(BAM_REMOVE_ORPHANS.out.versions.first())

    //
    // Sort, index PE BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS (
        BAM_REMOVE_ORPHANS.out.bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions.first())

    emit:
    name_bam = SAMTOOLS_SORT.out.bam                                                     // channel: [ val(meta), [ bam ] ]
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam.mix(ch_bam.single_end)                    // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai.mix(SAMTOOLS_INDEX.out.bai)               // channel: [ val(meta), [ bai ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats.mix(BAM_STATS_SAMTOOLS.out.stats)       // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat.mix(BAM_STATS_SAMTOOLS.out.flagstat) // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats.mix(BAM_STATS_SAMTOOLS.out.idxstats) // channel: [ val(meta), [ idxstats ] ]
    versions = ch_versions                                                               // channel: [ versions.yml ]
}
