/*
 * Filter BAM file
 */

include { FILTER_BAM         } from '../process/filter_bam'
include { REMOVE_BAM_ORPHANS } from '../process/remove_bam_orphans'
include { SORT_BAM           } from '../../nf-core/subworkflow/sort_bam'

workflow CLEAN_BAM {
    take:
    ch_bam_bai          // channel: [ val(meta), [ bam ], [bai] ]
    ch_bed              // channel: [ bed ]
    config              //    file: BAMtools filter JSON config file
    filter_bam_opts     //     map: options for filter_bam module
    rm_bam_orphans_opts //     map: options for remove_bam_orphans module
    samtools_opts       //     map: options for SAMTools modules

    main:
    FILTER_BAM(ch_bam_bai, ch_bed, config, filter_bam_opts)
    REMOVE_BAM_ORPHANS(FILTER_BAM.out.bam, rm_bam_orphans_opts)
    SORT_BAM(REMOVE_BAM_ORPHANS.out.bam, samtools_opts)

    emit:
    name_bam = REMOVE_BAM_ORPHANS.out.bam     // channel: [ val(meta), [ bam ] ]
    bam = SORT_BAM.out.bam                    // channel: [ val(meta), [ bam ] ]
    bai = SORT_BAM.out.bai                    // channel: [ val(meta), [ bai ] ]
    stats = SORT_BAM.out.stats                // channel: [ val(meta), [ stats ] ]
    flagstat = SORT_BAM.out.flagstat          // channel: [ val(meta), [ flagstat ] ]
    idxstats = SORT_BAM.out.idxstats          // channel: [ val(meta), [ idxstats ] ]
    bamtools_version = FILTER_BAM.out.version //    path: *.version.txt
}
