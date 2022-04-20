/*
 * Consensus peaks across samples, create boolean filtering file, SAF file for featureCounts
 */
process MACS2_CONSENSUS {
    tag "$meta.id"
    label 'process_long'

    conda (params.enable_conda ? "conda-forge::biopython conda-forge::r-optparse=1.7.1 conda-forge::r-upsetr=1.4.0 bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2f48cc59b03027e31ead6d383fe1b8057785dd24:5d182f583f4696f4c4d9f3be93052811b383341f-0':
        'quay.io/biocontainers/mulled-v2-2f48cc59b03027e31ead6d383fe1b8057785dd24:5d182f583f4696f4c4d9f3be93052811b383341f-0' }"

    input:
    tuple val(meta), path(peaks)

    output:
    tuple val(meta), path("*.bed")          , emit: bed
    tuple val(meta), path("*.saf")          , emit: saf
    tuple val(meta), path("*.pdf")          , emit: pdf
    tuple val(meta), path("*.boolean.txt")  , emit: boolean_txt
    tuple val(meta), path("*.intersect.txt"), emit: intersect_txt
    path "versions.yml"                     , emit: versions

    when:
    (meta.multiple_groups || meta.replicates_exist)

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/

    def prefix       = task.ext.prefix    ?: "${meta.id}"
    def peak_type    = params.narrow_peak ? 'narrowPeak' : 'broadPeak'
    def mergecols    = params.narrow_peak ? (2..10).join(',') : (2..9).join(',')
    def collapsecols = params.narrow_peak ? (['collapse']*9).join(',') : (['collapse']*8).join(',')
    def expandparam  = params.narrow_peak ? '--is_narrow_peak' : ''
    """
    sort -T '.' -k1,1 -k2,2n ${peaks.collect{it.toString()}.sort().join(' ')} \\
        | mergeBed -c $mergecols -o $collapsecols > ${prefix}.txt

    macs2_merged_expand.py \\
        ${prefix}.txt \\
        ${peaks.collect{it.toString()}.sort().join(',').replaceAll("_peaks.${peak_type}","")} \\
        ${prefix}.boolean.txt \\
        --min_replicates $params.min_reps_consensus \\
        $expandparam

    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2, \$3, \$4, "0", "+" }' ${prefix}.boolean.txt > ${prefix}.bed

    echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${prefix}.saf
    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$4, \$1, \$2, \$3,  "+" }' ${prefix}.boolean.txt >> ${prefix}.saf

    plot_peak_intersect.r -i ${prefix}.boolean.intersect.txt -o ${prefix}.boolean.intersect.plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """

}
