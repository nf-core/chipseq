/*
 * Get chromosome sizes from a fasta file
 */
 process GET_CHROM_SIZES {
     tag "$fasta"
     publishDir "${params.outdir}/genome", mode: params.publish_dir_mode

     input:
     path fasta

     output:
     path '*.sizes', emit: sizes
     path '*.fai', emit: fai

     script:
     """
     samtools faidx $fasta
     cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
     """
 }
