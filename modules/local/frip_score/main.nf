process FRIP_SCORE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92ad9212aec8a409dba3ed5468c74d59c93f4eb30a663379f8b405e24fbbf97f/data':
        'community.wave.seqera.io/library/bedtools_samtools:163d535533d93abd' }"

    input:
    tuple val(meta), path(bam), path(peak)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//'"), topic: versions, emit: versions_bedtools
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    READS_IN_PEAKS=\$(intersectBed -a $bam -b $peak $args | awk -F '\t' '{sum += \$NF} END {print sum}')
    samtools flagstat $bam > ${bam}.flagstat
    grep 'mapped (' ${bam}.flagstat | grep -v "primary" | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${prefix}", a/\$1}' > ${prefix}.FRiP.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.FRiP.txt
    """
}
