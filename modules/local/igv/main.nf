/*
 * Create IGV session file
 */
process IGV {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/62/622d8944750bc95bb56b4c3ed5c2b827e677c14073d48a5231e0f2bec0718add/data' :
        'community.wave.seqera.io/library/python:3.12.12--74abbf3898230efd' }"

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
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

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
    """
}
