/*
 * Map read(s) with bwa mem
 */
process BWA_MEM {
    tag "$name"
    label 'process_high'

    input:
    tuple val(name), val(single_end), path(reads)
    path index

    output:
    tuple val(name), val(single_end), path('*.bam')

    script:
    prefix = "${name}.Lb"
    rg = "\'@RG\\tID:${name}\\tSM:${name.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${name}\\tPU:1\'"
    if (params.seq_center) {
        rg = "\'@RG\\tID:${name}\\tSM:${name.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${name}\\tPU:1\\tCN:${params.seq_center}\'"
    }
    score = params.bwa_min_score ? "-T ${params.bwa_min_score}" : ''
    """
    bwa mem \\
        -t $task.cpus \\
        -M \\
        -R $rg \\
        $score \\
        ${index}/${params.bwa_base} \\
        $reads \\
        | samtools view -@ $task.cpus -b -h -F 0x0100 -O BAM -o ${prefix}.bam -
    """
}
