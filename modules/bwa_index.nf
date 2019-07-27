// Initialise parameters
params.outdir = './results'
params.saveGenomeIndex = false

process bwa_index {
    tag "$fasta"
    //label 'process_high'
    publishDir path: { params.saveGenomeIndex ? "${params.outdir}/reference_genome" : params.outdir },
               saveAs: { params.saveGenomeIndex ? it : null }, mode: 'copy'

    input:
    file fasta

    output:
    file "BWAIndex"

    script:
    """
    bwa index -a bwtsw $fasta
    mkdir BWAIndex && mv ${fasta}* BWAIndex
    """
}
