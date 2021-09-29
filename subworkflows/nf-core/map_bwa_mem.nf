/*
 * Map reads, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

params.bwa_mem_options  = [:]
params.samtools_sort_options  = [:]
params.samtools_index_options = [:]
params.samtools_stats_options = [:]

include { BWA_MEM           } from '../../modules/nf-core/modules/bwa/mem/main' addParams( options: params.bwa_mem_options )
include { BAM_SORT_SAMTOOLS } from './bam_sort_samtools'                        addParams( sort_options: params.samtools_sort_options, index_options: params.samtools_index_options, stats_options: params.samtools_stats_options )

workflow MAP_BWA_MEM {
    take:
    ch_reads         // channel: [ val(meta), [ reads ] ]
    ch_index         //    path: /path/to/index

    main:
    BWA_MEM(ch_reads, ch_index)
    BAM_SORT_SAMTOOLS(BWA_MEM.out.bam)

    emit:
    bam              = BAM_SORT_SAMTOOLS.out.bam              // channel: [ val(meta), [ bam ] ]
    bai              = BAM_SORT_SAMTOOLS.out.bai              // channel: [ val(meta), [ bai ] ]
    stats            = BAM_SORT_SAMTOOLS.out.stats            // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat         // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats         // channel: [ val(meta), [ idxstats ] ]
    bwa_version      = BWA_MEM.out.version                    //    path: versions.yml
    samtools_version = BAM_SORT_SAMTOOLS.out.samtools_version //    path: versions.yml
}
