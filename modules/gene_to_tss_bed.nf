// Initialise parameters
params.outdir = './results'

/*
 * PREPROCESSING - Generate TSS BED file
 */
process gene_to_tss_bed {
    tag "$bed"
    publishDir "${params.outdir}/reference_genome", mode: 'copy'

    input:
    file bed

    output:
    file "*.bed"

    script:
    """
    cat $bed | awk -v FS='\t' -v OFS='\t' '{ if(\$6=="+") \$3=\$2+1; else \$2=\$3-1; print \$1, \$2, \$3, \$4, \$5, \$6;}' > ${bed.baseName}.tss.bed
    """
}
