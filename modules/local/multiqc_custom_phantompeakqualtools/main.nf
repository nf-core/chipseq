process MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS {
    tag "$meta.id"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c73daa0b0040137fbea15fc8ee54c9bf4e0a9c5e9412cc7c13f7b38cc9c8bbd9/data':
        'community.wave.seqera.io/library/r-base:4.5.3--6814a4ccafc04d08' }"

    input:
    tuple val(meta), path(spp), path(rdata)
    path nsc_header
    path rsc_header
    path correlation_header

    output:
    tuple val(meta), path("*.spp_nsc_mqc.tsv")        , emit: nsc
    tuple val(meta), path("*.spp_rsc_mqc.tsv")        , emit: rsc
    tuple val(meta), path("*.spp_correlation_mqc.tsv"), emit: correlation
    tuple val("${task.process}"), val('R'), eval('R --version | sed "1!d; s/.*version //; s/ .*//"'), topic: versions, emit: versions_r

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp $correlation_header ${prefix}.spp_correlation_mqc.tsv
    Rscript --max-ppsize=500000 -e "load('$rdata'); write.table(crosscorr\\\$cross.correlation, file=\\"${prefix}.spp_correlation_mqc.tsv\\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"

    awk -v OFS='\t' '{print "${meta.id}", \$9}'  $spp | cat $nsc_header - > ${prefix}.spp_nsc_mqc.tsv
    awk -v OFS='\t' '{print "${meta.id}", \$10}' $spp | cat $rsc_header - > ${prefix}.spp_rsc_mqc.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.spp_nsc_mqc.tsv
    touch ${prefix}.spp_rsc_mqc.tsv
    touch ${prefix}.spp_correlation_mqc.tsv
    """
}
