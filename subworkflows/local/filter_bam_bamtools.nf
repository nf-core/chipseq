/*
 * Filter BAM file
 */

include { BAM_FILTER         } from '../../modules/local/bam_filter'
include { BAM_REMOVE_ORPHANS } from '../../modules/local/bam_remove_orphans'
include { BAM_SORT_SAMTOOLS  } from '../nf-core/bam_sort_samtools'

workflow FILTER_BAM_BAMTOOLS {
    take:
    ch_bam_bai                // channel: [ val(meta), [ bam ], [bai] ]
    ch_bed                    // channel: [ bed ]
    bamtools_filter_se_config //    file: BAMtools filter JSON config file for SE data
    bamtools_filter_pe_config //    file: BAMtools filter JSON config file for PE data

    main:
    ch_versions = Channel.empty()

    BAM_FILTER(ch_bam_bai, ch_bed, bamtools_filter_se_config, bamtools_filter_pe_config)
    BAM_REMOVE_ORPHANS(BAM_FILTER.out.bam)
    BAM_SORT_SAMTOOLS(BAM_REMOVE_ORPHANS.out.bam)

    ch_versions = ch_versions.mix(BAM_FILTER.out.versions,
                    BAM_REMOVE_ORPHANS.out.versions,
                    BAM_SORT_SAMTOOLS.out.versions)

    emit:
    name_bam = BAM_REMOVE_ORPHANS.out.bam     // channel: [ val(meta), [ bam ] ]
    bam      = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats    = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}
