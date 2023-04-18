/*
 * Map reads, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { STAR_ALIGN              } from '../../modules/local/star_align'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'

workflow ALIGN_STAR {
    take:
    ch_reads   // channel: [ val(meta), [ reads ] ]
    ch_index   // channel: /path/to/star/index/
    ch_fasta   // channel: /path/to/fasta
    seq_center //  string: sequencing center

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with STAR
    //
    STAR_ALIGN ( ch_reads, ch_index, seq_center )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( STAR_ALIGN.out.bam, ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    orig_bam       = STAR_ALIGN.out.bam                   // channel: [ val(meta), bam            ]
    log_final      = STAR_ALIGN.out.log_final             // channel: [ val(meta), log_final      ]
    log_out        = STAR_ALIGN.out.log_out               // channel: [ val(meta), log_out        ]
    log_progress   = STAR_ALIGN.out.log_progress          // channel: [ val(meta), log_progress   ]
    bam_sorted     = STAR_ALIGN.out.bam_sorted            // channel: [ val(meta), bam_sorted     ]
    bam_transcript = STAR_ALIGN.out.bam_transcript        // channel: [ val(meta), bam_transcript ]
    fastq          = STAR_ALIGN.out.fastq                 // channel: [ val(meta), fastq          ]
    tab            = STAR_ALIGN.out.tab                   // channel: [ val(meta), tab            ]

    bam            = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai            = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats          = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat       = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats       = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions       = ch_versions                          // channel: [ versions.yml ]
}
