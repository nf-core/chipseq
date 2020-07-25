process PHANTOMPEAKQUALTOOLS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else if (filename.endsWith('.version.txt')) null
                      else filename }

    container "quay.io/biocontainers/phantompeakqualtools:1.2.2--0"
    //container "https://depot.galaxyproject.org/singularity/phantompeakqualtools:1.2.2--0"

    conda (params.conda ? "bioconda::phantompeakqualtools=1.2.2" : null)

    input:
    tuple val(meta), path(bam)
    val opts

    output:
    tuple val(meta), path("*.out"), emit: spp
    tuple val(meta), path("*.pdf"), emit: pdf
    path "*.version.txt", emit: version

    script:
    prefix = opts.suffix ? "${meta.id}${opts.suffix}" : "${meta.id}"
    """
    RUN_SPP=`which run_spp.R`
    Rscript -e "library(caTools); source(\\"\$RUN_SPP\\")" -c="$bam" -savp="${prefix}.spp.pdf" -savd="${prefix}.spp.Rdata" -out="${prefix}.spp.out" -p=$task.cpus
    Rscript -e "load('${prefix}.spp.Rdata'); write.table(crosscorr\\\$cross.correlation, file=\\"${prefix}_spp_correlation_mqc.tsv\\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"

    touch phantompeakqualtools.version.txt
    """
}
