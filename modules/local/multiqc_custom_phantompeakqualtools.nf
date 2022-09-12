process MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS {
    conda (params.enable_conda ? "conda-forge::r-base=3.5.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:3.5.1':
        'quay.io/biocontainers/r-base:3.5.1' }"

    input:
    tuple val(meta), path(spp), path(rdata)
    path nsc_header
    path rsc_header
    path correlation_header

    output:
    tuple val(meta), path("*.spp_nsc_mqc.tsv")        , emit: nsc
    tuple val(meta), path("*.spp_rsc_mqc.tsv")        , emit: rsc
    tuple val(meta), path("*.spp_correlation_mqc.tsv"), emit: correlation

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp $correlation_header ${prefix}.spp_correlation_mqc.tsv
    Rscript --max-ppsize=500000 -e "load('$rdata'); write.table(crosscorr\\\$cross.correlation, file=\\"${prefix}.spp_correlation_mqc.tsv\\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"

    awk -v OFS='\t' '{print "${meta.id}", \$9}'  $spp | cat $nsc_header - > ${prefix}.spp_nsc_mqc.tsv
    awk -v OFS='\t' '{print "${meta.id}", \$10}' $spp | cat $rsc_header - > ${prefix}.spp_rsc_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
