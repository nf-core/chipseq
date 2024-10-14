/*
 * Create IGV session file
 */
process IGV {

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3':
        'biocontainers/python:3.8.3' }"

    input:
    val aligner_dir
    val peak_dir
    path fasta
    path ("${aligner_dir}/merged_library/bigwig/*")
    path ("${aligner_dir}/merged_library/macs3/${peak_dir}/*")
    path ("${aligner_dir}/merged_library/macs3/${peak_dir}/consensus/*")
    path ("mappings/*")

    output:
    path "*files.txt"  , emit: txt
    path "*.xml"       , emit: xml
    path fasta         , emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // scripts are bundled with the pipeline in nf-core/chipseq/bin/
    def consensus_dir = "${aligner_dir}/merged_library/macs3/${peak_dir}/consensus/*"
    """
    find * -type l -name "*.bigWig" -exec echo -e ""{}"\\t0,0,178" \\; > bigwig.igv.txt
    find * -type l -name "*Peak" -exec echo -e ""{}"\\t0,0,178" \\; > peaks.igv.txt
    # Avoid error when consensus not produced
    find * -type l -name "*.bed" -exec echo -e ""{}"\\t0,0,178" \\; | { grep "^$consensus_dir" || test \$? = 1; } > consensus.igv.txt

    if [ -d "mappings" ]; then
        cat mappings/* > replace_paths.txt
    else
        touch replace_paths.txt
    fi

    cat *.igv.txt > igv_files_orig.txt
    igv_files_to_session.py igv_session.xml igv_files_orig.txt replace_paths.txt ../../genome/${fasta.getName()} --path_prefix '../../'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
