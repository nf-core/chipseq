process DESEQ2_QC {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::r-base=4.0 bioconda::bioconductor-deseq2=1.28.0 bioconda::bioconductor-biocparallel bioconda::bioconductor-tximport bioconda::bioconductor-complexheatmap conda-forge::r-optparse conda-forge::r-ggplot2 conda-forge::r-rcolorbrewer conda-forge::r-pheatmap" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
    tuple val(meta), path(counts)
    path deseq2_pca_header
    path deseq2_clustering_header

    output:
    path "*.pdf"                , optional:true, emit: pdf
    path "*.RData"              , optional:true, emit: rdata
    path "*pca.vals.txt"        , optional:true, emit: pca_txt
    path "*pca.vals_mqc.tsv"    , optional:true, emit: pca_multiqc
    path "*sample.dists.txt"    , optional:true, emit: dists_txt
    path "*sample.dists_mqc.tsv", optional:true, emit: dists_multiqc
    path "*.log"                , optional:true, emit: log
    path "size_factors"         , optional:true, emit: size_factors
    path "versions.yml"         , emit: versions

    script:
    def args      = task.ext.args ?: ''
    def peak_type = params.narrow_peak ? 'narrowPeak' : 'broadPeak'
    def antibody  = meta.antibody
    def prefix    = "${antibody}.consensus_peaks"
    """
    deseq2_qc.r \\
        --count_file $counts \\
        --outdir ./ \\
        --outprefix $prefix \\
        --cores $task.cpus \\
        $args

    sed 's/deseq2_pca/deseq2_pca_${task.index}/g' <$deseq2_pca_header >tmp.txt
    sed -i -e 's/DESeq2 /${antibody} DESeq2 /g' tmp.txt
    cat tmp.txt ${prefix}.pca.vals.txt > ${prefix}.pca.vals_mqc.tsv

    sed 's/deseq2_clustering/deseq2_clustering_${task.index}/g' <$deseq2_clustering_header >tmp.txt
    sed -i -e 's/DESeq2 /${antibody} DESeq2 /g' tmp.txt
    cat tmp.txt ${prefix}.sample.dists.txt > ${prefix}.sample.dists_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
