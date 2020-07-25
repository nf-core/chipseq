/*
 * Map reads, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { BWA_MEM  } from '../software/bwa_mem'
include { SORT_BAM } from './sort_bam'

workflow MAP_READS {
    take:
    ch_reads      // channel: [ val(meta), [ reads ] ]
    ch_index      //    path: /path/to/index
    ch_fasta      //    path: /path/to/genome.fasta
    bwa_mem_opts  //     map: options for BWA MEM module
    samtools_opts //     map: options for SAMTools modules

    main:
    BWA_MEM(ch_reads, ch_index, ch_fasta, bwa_mem_opts)
    SORT_BAM(BWA_MEM.out.bam, samtools_opts)

    emit:
    bam = SORT_BAM.out.bam                           // channel: [ val(meta), [ bam ] ]
    bai = SORT_BAM.out.bai                           // channel: [ val(meta), [ bai ] ]
    stats = SORT_BAM.out.stats                       // channel: [ val(meta), [ stats ] ]
    flagstat = SORT_BAM.out.flagstat                 // channel: [ val(meta), [ flagstat ] ]
    idxstats = SORT_BAM.out.idxstats                 // channel: [ val(meta), [ idxstats ] ]
    bwa_version = BWA_MEM.out.version                //    path: *.version.txt
    samtools_version = SORT_BAM.out.samtools_version //    path: *.version.txt
}
