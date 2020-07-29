def SOFTWARE = 'bwa'

process BWA_MEM {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/${options.publish_dir}${options.publish_by_id ? "/${meta.id}" : ''}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (options.publish_results == "none") null
                      else if (filename.endsWith('.version.txt')) null
                      else filename }

    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:eabfac3657eda5818bae4090db989e3d41b01542-0"
    //container "https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:eabfac3657eda5818bae4090db989e3d41b01542-0"

    conda (params.conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.10" : null)

    input:
    tuple val(meta), path(reads)
    path index
    path fasta
    val options

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.version.txt", emit: version

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    rg = meta.read_group ? "-R ${meta.read_group}" : ""
    """
    bwa mem \\
        $options.args \\
        $rg \\
        -t $task.cpus \\
        $fasta \\
        $reads \\
        | samtools view $options.args2 -@ $task.cpus -bS -o ${prefix}.bam -

    echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//' > ${SOFTWARE}.version.txt
    """
}
