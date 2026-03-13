/*
 * Prepare genome intervals for filtering by removing regions in blacklist file
 */
process GENOME_BLACKLIST_REGIONS {
    tag "$sizes"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/6397750e9730a3fbcc5b4c43f14bd141c64c723fd7dad80e47921a68a7c3cd21/data'
        : 'community.wave.seqera.io/library/bedtools_coreutils:a623c13f66d5262b'}"

    input:
    path sizes
    path blacklist

    output:
    path '*.bed'       , emit: bed
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def file_out = "${sizes.simpleName}.include_regions.bed"
    if (blacklist) {
        """
        sortBed -i $blacklist -g $sizes | complementBed -i stdin -g $sizes > $file_out
        """
    } else {
        """
        awk '{print \$1, '0' , \$2}' OFS='\t' $sizes > $file_out
        """
    }
}
