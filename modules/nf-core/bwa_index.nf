/*
 * Build BWA index
 */
process BWA_INDEX {
    tag "$fasta"
    //label 'process_high'
    // publishDir path: { params.save_reference ? "${params.outdir}/genome" : params.outdir },
    //     saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

    input:
    path fasta

    output:
    path 'BWAIndex'

    script:
    """
    bwa index -a bwtsw $fasta
    mkdir BWAIndex && mv ${fasta}* BWAIndex
    """
}
