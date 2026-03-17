/*
 * Consensus peaks across samples, create boolean filtering file, SAF file for featureCounts
 */
process MACS3_CONSENSUS {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/57/573950f2d84f1e424be5344ed579deb5082fd23751e4a6142d7c584bb3af2d3e/data':
        'community.wave.seqera.io/library/bedtools_biopython_r-optparse_r-upsetr:c682828913318fba' }"

    input:
    tuple val(meta), path(peaks)
    val is_narrow_peak

    output:
    tuple val(meta), path("*.bed")          , emit: bed
    tuple val(meta), path("*.saf")          , emit: saf
    tuple val(meta), path("*.pdf")          , emit: pdf
    tuple val(meta), path("*.antibody.txt") , emit: txt
    tuple val(meta), path("*.boolean.txt")  , emit: boolean_txt
    tuple val(meta), path("*.intersect.txt"), emit: intersect_txt

    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('R'), eval("R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//'"), topic: versions, emit: versions_r

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def peak_type    = is_narrow_peak  ? 'narrowPeak' : 'broadPeak'
    def mergecols    = is_narrow_peak  ? (2..10).join(',') : (2..9).join(',')
    def collapsecols = is_narrow_peak  ? (['collapse']*9).join(',') : (['collapse']*8).join(',')
    def expandparam  = is_narrow_peak  ? '--is_narrow_peak' : ''
    """
    sort -T '.' -k1,1 -k2,2n ${peaks.collect{it.toString()}.sort().join(' ')} \\
        | mergeBed -c $mergecols -o $collapsecols > ${prefix}.txt

    macs3_merged_expand.py \\
        ${prefix}.txt \\
        ${peaks.collect{it.toString()}.sort().join(',').replaceAll("_peaks.${peak_type}","")} \\
        ${prefix}.boolean.txt \\
        --min_replicates $params.min_reps_consensus \\
        $args \\
        $expandparam

    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2, \$3, \$4, "0", "+" }' ${prefix}.boolean.txt > ${prefix}.bed

    echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${prefix}.saf
    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$4, \$1, \$2, \$3,  "+" }' ${prefix}.boolean.txt >> ${prefix}.saf

    plot_peak_intersect.r -i ${prefix}.boolean.intersect.txt -o ${prefix}.boolean.intersect.plot.pdf

    echo "${prefix}.bed\t${meta.id}/${prefix}.bed" > ${prefix}.antibody.txt
    """
    stub:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.antibody.txt
    """
}
