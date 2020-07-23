process BWA_MEM {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else if (filename.endsWith('.version.txt')) null
                    else filename }

    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:eabfac3657eda5818bae4090db989e3d41b01542-0"
    //container "https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:eabfac3657eda5818bae4090db989e3d41b01542-0"

    conda (params.conda ? "${moduleDir}/environment.yml" : null)

    input:
    tuple val(meta), path(reads)
    path index
    path fasta
    val opts

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.version.txt", emit: version

    script:
    prefix = opts.suffix ? "${meta.id}${opts.suffix}" : "${meta.id}"
    rg = meta.read_group ? "-R ${meta.read_group}" : ""
    """
    bwa mem \\
        $opts.args \\
        $rg \\
        -t $task.cpus \\
        $fasta \\
        $reads \\
        | samtools view $opts.args2 -@ $task.cpus -bS -o ${prefix}.bam -
    echo \$(bwa 2>&1) | sed -n "s/.*\\(v.*\$\\)/\\1/p" > bwa.version.txt
    """
}
