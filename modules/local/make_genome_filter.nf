// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]

/*
 * Prepare genome intervals for filtering by removing regions in blacklist file
 */
process MAKE_GENOME_FILTER {
    tag "$sizes"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'genome', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0"
    } else {
        container "quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"
    }

    input:
    path sizes
    path blacklist

    output:
    path '*.bed'       , emit: bed
    path "versions.yml", emit: version

    script:
    def software = 'bedtools'
    def file_out = "${sizes.simpleName}.include_regions.bed"
    if (params.blacklist) {
        """
        sortBed -i $blacklist -g $sizes | complementBed -i stdin -g $sizes > $file_out

        cat <<-END_VERSIONS > versions.yml
            ${getProcessName(task.process)}:
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    } else {
        """
        awk '{print \$1, '0' , \$2}' OFS='\t' $sizes > $file_out

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    }
}
