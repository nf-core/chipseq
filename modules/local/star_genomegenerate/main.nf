process STAR_GENOMEGENERATE {
    tag "$fasta"
    label 'process_high'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52097ce908b427400215cbcfbf143810dc2e901b043cc153ea1d4ef5b2b025e9/data' :
        'community.wave.seqera.io/library/samtools_star_gawk:da63acc62f0ad3b8' }"

    input:
    path fasta
    path gtf

    output:
    path "star"        , emit: index
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools
    tuple val("${task.process}"), val('star'), eval('STAR --version | sed -e "s/STAR_//"'), topic: versions, emit: versions_star

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = (task.ext.args ?: '').tokenize()
    def memory = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    if (args.contains('--genomeSAindexNbases')) {
        """
        mkdir star
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta \\
            --sjdbGTFfile $gtf \\
            --runThreadN $task.cpus \\
            $memory \\
            ${args.join(' ')}
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
        END_VERSIONS
        """
    } else {
        """
        samtools faidx $fasta
        NUM_BASES=`gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${fasta}.fai`
        mkdir star
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta \\
            --sjdbGTFfile $gtf \\
            --runThreadN $task.cpus \\
            --genomeSAindexNbases \$NUM_BASES \\
            $memory \\
            ${args.join(' ')}
        """
    }
}
