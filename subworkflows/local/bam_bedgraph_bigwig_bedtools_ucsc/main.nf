//
// Convert BAM to normalised bigWig via bedGraph using BEDTools and UCSC
//

include { BEDTOOLS_GENOMECOV    } from '../../../modules/nf-core/bedtools/genomecov'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/nf-core/ucsc/bedgraphtobigwig'

workflow BAM_BEDGRAPH_BIGWIG_BEDTOOLS_UCSC {
    take:
    ch_bam_flagstat // channel: [ val(meta), [bam], [flagstat] ]
    ch_chrom_sizes  // channel: [ bed ]

    main:


    //
    // Extract scale factor from flagstat and prepare input for bedtools genomecov
    //
    ch_bam_scale = ch_bam_flagstat
        .map { meta, bam, flagstat ->
            // Parse flagstat to get mapped reads count
            def flagstat_content = flagstat.text
            def match = flagstat_content =~ /(\d+) \+ \d+ mapped/
            if (match.size() == 0) {
                // Try alternative pattern for different flagstat formats or stub files
                match = flagstat_content =~ /(\d+) mapped/
                if (match.size() == 0) {
                    // For stub tests, use a default scale factor
                    log.warn "Could not parse mapped reads from flagstat, using default scale factor of 1.0"
                    return [meta, bam, 1.0]
                }
            }
            def mapped_reads = match[0][1] as Integer
            if (mapped_reads == 0) {
                log.warn "Zero mapped reads found in flagstat file, using default scale factor of 1.0"
                return [meta, bam, 1.0]
            }
            def scale_factor = 1000000 / mapped_reads
            [meta, bam, scale_factor]
        }

    //
    // Create bedGraph coverage track
    //
    BEDTOOLS_GENOMECOV (
        ch_bam_scale,
        [],
        'bedGraph',
        true
    )

    //
    // Create bigWig coverage tracks
    //
    UCSC_BEDGRAPHTOBIGWIG (
        BEDTOOLS_GENOMECOV.out.genomecov,
        ch_chrom_sizes
    )

    emit:
    bedgraph     = BEDTOOLS_GENOMECOV.out.genomecov     // channel: [ val(meta), [ bedgraph ] ]
    bigwig       = UCSC_BEDGRAPHTOBIGWIG.out.bigwig     // channel: [ val(meta), [ bigwig ] ]

}
