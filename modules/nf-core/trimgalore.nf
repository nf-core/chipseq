/*
 * Trim Galore!
 */
process TRIMGALORE {
    tag "$name"
    label 'process_high'
    publishDir "${params.outdir}/trim_galore", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith('.html')) "fastqc/$filename"
                      else if (filename.endsWith('.zip')) "fastqc/zips/$filename"
                      else if (filename.endsWith('trimming_report.txt')) "logs/$filename"
                      else params.save_trimmed ? filename : null
                }

    input:
    tuple val(name), val(single_end), path(reads)

    output:
    tuple val(name), val(single_end), path('*.fq.gz'), emit: reads
    path '*.txt', emit: log
    path '*.{zip,html}', emit: fastqc

    script:
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (single_end) cores = (task.cpus as int) - 3
        if (cores < 1) cores = 1
        if (cores > 4) cores = 4
    }

    c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
    nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''

    // Added soft-links to original fastqs for consistent naming in MultiQC
    if (single_end) {
        """
        [ ! -f  ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
        trim_galore --cores $cores --fastqc --gzip $c_r1 $tpc_r1 $nextseq ${name}.fastq.gz
        """
    } else {
        """
        [ ! -f  ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
        [ ! -f  ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
        trim_galore --cores $cores --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $nextseq ${name}_1.fastq.gz ${name}_2.fastq.gz
        """
    }
}
