// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    //TODO substitute with a newest tag (see https://github.com/BioContainers/containers/issues/416)
    conda (params.enable_conda ? "conda-forge::r-base=3.5.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/r-base:3.5.1"
    } else {
        container "quay.io/biocontainers/r-base:3.5.1"
    }

    input:
    tuple val(meta), path(spp), path(rdata)
    path nsc_header
    path rsc_header
    path correlation_header

    output:
    tuple val(meta), path("*.spp_nsc_mqc.tsv"), emit: nsc
    tuple val(meta), path("*.spp_rsc_mqc.tsv"), emit: rsc
    tuple val(meta), path("*.spp_correlation_mqc.tsv"), emit: correlation

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    cp $correlation_header ${prefix}.spp_correlation_mqc.tsv
    Rscript -e "load('$rdata'); write.table(crosscorr\\\$cross.correlation, file=\\"${prefix}.spp_correlation_mqc.tsv\\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"

    awk -v OFS='\t' '{print "${meta.id}", \$9}'  $spp | cat $nsc_header - > ${prefix}.spp_nsc_mqc.tsv
    awk -v OFS='\t' '{print "${meta.id}", \$10}' $spp | cat $rsc_header - > ${prefix}.spp_rsc_mqc.tsv
    """
}
