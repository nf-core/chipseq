include { SAMTOOLS_SORT           } from '../../../modules/nf-core/samtools/sort'
include { SAMTOOLS_INDEX          } from '../../../modules/nf-core/samtools/index'
include { BAM_SORT_STATS_SAMTOOLS } from '../../nf-core/bam_sort_stats_samtools'
include { BAM_STATS_SAMTOOLS      } from '../../nf-core/bam_stats_samtools'

include { BAMTOOLS_FILTER         } from '../../../modules/local/bamtools_filter'
include { BAM_REMOVE_ORPHANS      } from '../../../modules/local/bam_remove_orphans'

workflow BAM_FILTER_BAMTOOLS {
    take:
    ch_bam_bai                   // channel: [ val(meta), [ bam ], [bai] ]
    ch_bed                       // channel: [ bed ]
    ch_fasta                     // channel: [ val(meta), fasta ]
    ch_bamtools_filter_se_config // channel: [ config_file ]
    ch_bamtools_filter_pe_config // channel: [ config_file ]

    main:


    //
    // Filter BAM file with BAMTools
    //
    BAMTOOLS_FILTER (
        ch_bam_bai,
        ch_bed,
        ch_bamtools_filter_se_config,
        ch_bamtools_filter_pe_config
    )

    ch_bam_single_end = BAMTOOLS_FILTER.out.bam
        .filter { meta, _bam -> meta.single_end == true }
        .map { meta, bam -> [meta, bam] }

    ch_bam_paired_end = BAMTOOLS_FILTER.out.bam
        .filter { meta, _bam -> meta.single_end != true }
        .map { meta, bam -> [meta, bam] }

    //
    // Index SE BAM file
    //
    SAMTOOLS_INDEX (
        ch_bam_single_end
    )

    //
    // Run samtools stats, flagstat and idxstats on SE BAM
    //
    BAM_STATS_SAMTOOLS (
        ch_bam_single_end.join(SAMTOOLS_INDEX.out.index),
        ch_fasta.map { meta, fasta -> [ meta, fasta, [] ] }
    )

    //
    // Name sort PE BAM before filtering with pysam
    //
    SAMTOOLS_SORT (
        ch_bam_paired_end,
        ch_fasta.map { meta, fasta -> [ meta, fasta, [] ] },
        ''
    )

    //
    // Remove orphan reads from PE BAM file
    //
    BAM_REMOVE_ORPHANS (
        SAMTOOLS_SORT.out.bam
    )

    //
    // Sort, index PE BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS (
        BAM_REMOVE_ORPHANS.out.bam,
        ch_fasta.map { meta, fasta -> [ meta, fasta, [] ] }
    )

    emit:
    name_bam = SAMTOOLS_SORT.out.bam                                                       // channel: [ val(meta), [ bam ] ]
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam.mix(ch_bam_single_end)                      // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.index.mix(SAMTOOLS_INDEX.out.index)             // channel: [ val(meta), [ bai ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats.mix(BAM_STATS_SAMTOOLS.out.stats)         // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat.mix(BAM_STATS_SAMTOOLS.out.flagstat)   // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats.mix(BAM_STATS_SAMTOOLS.out.idxstats)   // channel: [ val(meta), [ idxstats ] ]
}
