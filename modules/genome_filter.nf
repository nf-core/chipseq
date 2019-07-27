// Initialise parameters
params.outdir = './results'
params.blacklist = false

process genome_filter {
    tag "$fasta"
    publishDir "${params.outdir}/reference_genome", mode: 'copy'

    input:
    file fasta

    output:
    file "$fasta"       // FASTA FILE FOR IGV
    file "*.fai"        // FAI INDEX FOR REFERENCE GENOME
    file "*.bed"        // BED FILE WITHOUT BLACKLIST REGIONS
    file "*.sizes"      // CHROMOSOME SIZES FILE FOR BEDTOOLS

    script:
    blacklist_filter = params.blacklist ? "sortBed -i ${params.blacklist} -g ${fasta}.sizes | complementBed -i stdin -g ${fasta}.sizes" : "awk '{print \$1, '0' , \$2}' OFS='\t' ${fasta}.sizes"
    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
    $blacklist_filter > ${fasta}.include_regions.bed
    """
}
