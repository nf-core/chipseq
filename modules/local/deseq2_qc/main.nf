process DESEQ2_QC {
    tag "$meta.id"
    label 'process_medium'

    // (Bio)conda packages have intentionally not been pinned to a specific version
    // This was to avoid the pipeline failing due to package conflicts whilst creating the environment when using -profile conda
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
    tuple val(meta), path(counts)
    path deseq2_pca_header
    path deseq2_clustering_header

    output:
    path "*.pdf"                , optional:true, emit: pdf
    path "*.RData"              , optional:true, emit: rdata
    path "*.rds"                , optional:true, emit: rds
    path "*pca.vals.txt"        , optional:true, emit: pca_txt
    path "*pca.vals_mqc.tsv"    , optional:true, emit: pca_multiqc
    path "*sample.dists.txt"    , optional:true, emit: dists_txt
    path "*sample.dists_mqc.tsv", optional:true, emit: dists_multiqc
    path "*.log"                , optional:true, emit: log
    path "size_factors"         , optional:true, emit: size_factors
    tuple val("${task.process}"), val('R'), eval('R --version | sed "1!d; s/.*version //; s/ .*//"'), topic: versions, emit: versions_r
    tuple val("${task.process}"), val('DESeq2'), eval('Rscript -e "library(DESeq2); cat(as.character(packageVersion(\'DESeq2\')))"'), topic: versions, emit: versions_deseq2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    deseq2_qc.r \\
        --count_file $counts \\
        --outdir ./ \\
        --outprefix $prefix \\
        --cores $task.cpus \\
        $args

    sed 's/deseq2_pca/deseq2_pca_${task.index}/g' <$deseq2_pca_header >tmp.txt
    sed -i -e 's/DESeq2 /${meta.id} DESeq2 /g' tmp.txt
    cat tmp.txt ${prefix}.pca.vals.txt > ${prefix}.pca.vals_mqc.tsv

    set +C

    sed 's/deseq2_clustering/deseq2_clustering_${task.index}/g' <$deseq2_clustering_header >tmp.txt
    sed -i -e 's/DESeq2 /${meta.id} DESeq2 /g' tmp.txt
    cat tmp.txt ${prefix}.sample.dists.txt > ${prefix}.sample.dists_mqc.tsv
    """
}
