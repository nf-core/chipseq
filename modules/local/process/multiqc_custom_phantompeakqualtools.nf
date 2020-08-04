// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename=filename, options=options, publish_dir=getSoftwareName(task.process), publish_id=meta.id) }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    tuple val(meta), path(spp), path(rdata)
    path nsc_header
    path rsc_header
    path correlation_header
    val options

    output:
    tuple val(meta), path("*.spp_nsc_mqc.tsv"), emit: nsc
    tuple val(meta), path("*.spp_rsc_mqc.tsv"), emit: rsc
    tuple val(meta), path("*.spp_correlation_mqc.tsv"), emit: correlation

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    prefix = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    """
    cp $correlation_header ${prefix}.spp_correlation_mqc.tsv
    Rscript -e "load('$rdata'); write.table(crosscorr\\\$cross.correlation, file=\\"${prefix}.spp_correlation_mqc.tsv\\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"

    awk -v OFS='\t' '{print "${meta.id}", \$9}'  $spp | cat $nsc_header - > ${prefix}.spp_nsc_mqc.tsv
    awk -v OFS='\t' '{print "${meta.id}", \$10}' $spp | cat $rsc_header - > ${prefix}.spp_rsc_mqc.tsv
    """
}
