// Initialise parameters
params.outdir = './results'

process gtf_to_bed {
    tag "$gtf"
    //label 'process_low'
    publishDir "${params.outdir}/reference_genome", mode: 'copy'

    input:
    file gtf

    output:
    file "*.bed"

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    """
    gtf2bed $gtf > ${gtf.baseName}.bed
    """
}
