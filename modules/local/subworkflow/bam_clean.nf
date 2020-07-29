/*
 * Filter BAM file
 */

include { BAM_FILTER         } from '../process/bam_filter'
include { BAM_REMOVE_ORPHANS } from '../process/bam_remove_orphans'
include { BAM_SORT_SAMTOOLS  } from '../../nf-core/subworkflow/bam_sort_samtools'

workflow BAM_CLEAN {
    take:
    ch_bam_bai                 // channel: [ val(meta), [ bam ], [bai] ]
    ch_bed                     // channel: [ bed ]
    bamtools_filter_se_config  //    file: BAMtools filter JSON config file for SE data
    bamtools_filter_pe_config  //    file: BAMtools filter JSON config file for PE data
    bam_filter_options         //     map: options for bam_filter module
    bam_remove_orphans_options //     map: options for bam_remove_orphans module
    samtools_options           //     map: options for SAMTools modules

    main:
    BAM_FILTER(ch_bam_bai, ch_bed, bamtools_filter_se_config, bamtools_filter_pe_config, bam_filter_options)
    BAM_REMOVE_ORPHANS(BAM_FILTER.out.bam, bam_remove_orphans_options)
    BAM_SORT_SAMTOOLS(BAM_REMOVE_ORPHANS.out.bam, samtools_options)

    emit:
    name_bam = BAM_REMOVE_ORPHANS.out.bam     // channel: [ val(meta), [ bam ] ]
    bam = BAM_SORT_SAMTOOLS.out.bam           // channel: [ val(meta), [ bam ] ]
    bai = BAM_SORT_SAMTOOLS.out.bai           // channel: [ val(meta), [ bai ] ]
    stats = BAM_SORT_SAMTOOLS.out.stats       // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    bamtools_version = BAM_FILTER.out.version //    path: *.version.txt
}
