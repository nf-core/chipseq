/*
 * Prepare genome intervals for filtering
 * by removing regions in blacklist file
 */
process MAKE_GENOME_FILTER {
    tag "$sizes"
    publishDir "${params.outdir}/genome", mode: params.publish_dir_mode

    input:
    path sizes
    path blacklist

    output:
    path '*.bed'

    script:
    file_out = "${sizes.simpleName}.include_regions.bed"
    if (params.blacklist) {
        """
        sortBed -i $blacklist -g $sizes | complementBed -i stdin -g $sizes > $file_out
        """
    } else {
        """
        awk '{print \$1, '0' , \$2}' OFS='\t' $sizes > $file_out
        """
    }
}
