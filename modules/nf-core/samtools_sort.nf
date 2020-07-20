def MODULE = "samtools_sort"
params.publish_dir = MODULE
params.publish_results = "default"
params.suffix = ''
//params.samtools_sort_args = ''

process SAMTOOLS_SORT {
    publishDir "${params.outdir}/${params.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (params.publish_results == "none") null
                    else filename }

    container "docker.pkg.github.com/nf-core/$MODULE"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(name), val(single_end), path(bam)
    //val (samtools_sort_args)

    output:
    tuple val(name), val(single_end), path('*.bam')
    path "*.version.txt", emit: version

    script:
    prefix = "${name}.${suffix}"
    """
    samtools sort -@ $task.cpus -o ${prefix}.sorted.bam -T $name $bam
    samtools --version | sed -n "s/.*\\(v.*\$\\)/\\1/p" > samtools.version.txt
    """
}
