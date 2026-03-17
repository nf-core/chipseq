process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9b/9b8ecb2f9a77b5e7573ef6fae2f4c2e771064f7a129ed1329913c1025c33f365/data' :
        'community.wave.seqera.io/library/star:2.7.11b--822039d47adf19a7' }"

    input:
    tuple val(meta) , path(reads)
    path  index
    val seq_center

    output:
    tuple val(meta), path('*d.out.bam')       , emit: bam
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    tuple val("${task.process}"), val('star'), eval('STAR --version | sed -e "s/STAR_//"'), topic: versions, emit: versions_star

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')               , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seq_center_tag = seq_center ? "--outSAMattrRGline ID:$prefix 'CN:$seq_center' 'SM:$prefix'" : "--outSAMattrRGline ID:$prefix 'SM:$prefix'"
    def out_sam_type = (args.contains('--outSAMtype')) ? '' : '--outSAMtype BAM Unsorted'
    def mv_unsorted_bam = (args.contains('--outSAMtype BAM Unsorted SortedByCoordinate')) ? "mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam" : ''
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn $reads  \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        $out_sam_type \\
        $seq_center_tag \\
        $args
    $mv_unsorted_bam
    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi
    """
}
