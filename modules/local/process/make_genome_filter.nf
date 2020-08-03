// Import generic module functions
include { initOptions; saveFiles } from './functions'

/*
 * Prepare genome intervals for filtering by removing regions in blacklist file
 */
process MAKE_GENOME_FILTER {
    tag "$sizes"
    publishDir "${params.outdir}/${options.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename, options, "genome") }

    container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    //container "https://depot.galaxyproject.org/singularity/bedtools:2.29.2--hc088bd4_0"

    conda (params.conda ? "bioconda::bedtools=2.29.2" : null)

    input:
    path sizes
    path blacklist
    val options

    output:
    path '*.bed', emit: bed
    path "*.version.txt", emit: version

    script:
    def software = 'bedtools'
    file_out = "${sizes.simpleName}.include_regions.bed"
    if (params.blacklist) {
        """
        sortBed -i $blacklist -g $sizes | complementBed -i stdin -g $sizes > $file_out
        bedtools --version | sed -e "s/bedtools v//g" > ${software}.version.txt
        """
    } else {
        """
        awk '{print \$1, '0' , \$2}' OFS='\t' $sizes > $file_out
        bedtools --version | sed -e "s/bedtools v//g" > ${software}.version.txt
        """
    }
}
