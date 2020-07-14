/*
 * Filter BAM file
 */
process FILTER_BAM {
    tag "$name"
    label 'process_medium'
    publishDir path: "${params.outdir}/bwa/mergedLibrary", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (single_end || params.save_align_intermeds) {
                          if (filename.endsWith('.flagstat')) "samtools_stats/$filename"
                          else if (filename.endsWith('.idxstats')) "samtools_stats/$filename"
                          else if (filename.endsWith('.stats')) "samtools_stats/$filename"
                          else if (filename.endsWith('.sorted.bam')) filename
                          else if (filename.endsWith('.sorted.bam.bai')) filename
                          else null
                      }
                }

    input:
    tuple val(name), val(single_end), path(bam)
    path bed
    path filter_config

    output:
    tuple val(name), val(single_end), path('*.bam'), emit: bam
    tuple val(name), val(single_end), path('*.bai'), emit: bai
    tuple val(name), val(single_end), path('*.stats'), emit: stats
    tuple val(name), val(single_end), path('*.flagstat'), emit: flagstat
    tuple val(name), val(single_end), path('*.idxstats'), emit: idxstats

    script:
    prefix = single_end ? "${name}.mLb.clN" : "${name}.mLb.flT"
    filter_params = single_end ? '-F 0x004' : '-F 0x004 -F 0x0008 -f 0x001'
    dup_params = params.keep_dups ? '' : '-F 0x0400'
    multimap_params = params.keep_multi_map ? '' : '-q 1'
    blacklist_params = params.blacklist ? "-L $bed" : ''
    name_sort_bam = single_end ? '' : "samtools sort -n -@ $task.cpus -o ${prefix}.bam -T $prefix ${prefix}.sorted.bam"
    """
    samtools view \\
        $filter_params \\
        $dup_params \\
        $multimap_params \\
        $blacklist_params \\
        -b ${bam[0]} \\
        | bamtools filter \\
            -out ${prefix}.sorted.bam \\
            -script $filter_config
    samtools index ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
    samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
    samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats

    $name_sort_bam
    """
}
