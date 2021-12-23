process PLOT_MACS2_QC {
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::r-base=4.0.3 conda-forge::r-reshape2=1.4.4 conda-forge::r-optparse=1.6.6 conda-forge::r-ggplot2=3.3.3 conda-forge::r-scales=1.1.1 conda-forge::r-viridis=0.5.1 conda-forge::r-tidyverse=1.3.0 bioconda::bioconductor-biostrings=2.58.0 bioconda::bioconductor-complexheatmap=2.6.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0fb10cdfa7:05678da05b8e5a7a5130e90a9f9a6c585b965afa-0':
        'quay.io/biocontainers/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0fb10cdfa7:05678da05b8e5a7a5130e90a9f9a6c585b965afa-0' }"

    input:
    path peaks

    output:
    path '*.txt'       , emit: txt
    path '*.pdf'       , emit: pdf
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args      = task.ext.args ?: ''
    def peak_type = params.narrow_peak ? 'narrowPeak' : 'broadPeak'
    """
    plot_macs2_qc.r \\
        -i ${peaks.join(',')} \\
        -s ${peaks.join(',').replaceAll("_peaks.${peak_type}","")} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
