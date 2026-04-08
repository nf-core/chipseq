/*
 * Create IGV session file (XML for IGV Desktop, JSON for IGV.js / Seqera Data Explorer)
 */
process IGV {

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3':
        'biocontainers/python:3.8.3' }"

    input:
    val aligner_dir
    val peak_dir
    val genome_id
    path fasta
    path ("${aligner_dir}/merged_library/bigwig/*")
    path ("${aligner_dir}/merged_library/macs3/${peak_dir}/*")
    path ("${aligner_dir}/merged_library/macs3/${peak_dir}/consensus/*")
    path ("mappings/*")
    path ip_bams
    path ip_bais
    path control_bams
    path control_bais
    val  ip_sample_ids
    val  control_sample_ids

    output:
    path "*files.txt"       , emit: txt
    path "*.xml"            , emit: xml
    path "igv_session.json" , emit: json
    path fasta              , emit: fasta
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // scripts are bundled with the pipeline in nf-core/chipseq/bin/
    def consensus_dir = "${aligner_dir}/merged_library/macs3/${peak_dir}/consensus/*"

    // Build space-separated file lists for the Python script
    def ip_bam_args     = ip_bams.name     != 'NO_IP_BAMS'     ? "--bam_files ${ip_bams.collect{ it.name }.join(' ')}"         : ''
    def ip_bai_args     = ip_bais.name     != 'NO_IP_BAIS'     ? "--bai_files ${ip_bais.collect{ it.name }.join(' ')}"         : ''
    def ctrl_bam_args   = control_bams.name != 'NO_CTRL_BAMS'  ? "--control_bam_files ${control_bams.collect{ it.name }.join(' ')}" : ''
    def ctrl_bai_args   = control_bais.name != 'NO_CTRL_BAIS'  ? "--control_bai_files ${control_bais.collect{ it.name }.join(' ')}" : ''
    def sample_id_args  = ip_sample_ids     ? "--sample_ids ${ip_sample_ids.join(' ')}"      : ''
    def control_id_args = control_sample_ids ? "--control_ids ${control_sample_ids.join(' ')}" : ''
    """
    find * -type l -name "*.bigWig" -exec echo -e ""{}"\\t0,0,178" \\; > bigwig.igv.txt
    find * -type l -name "*Peak" -exec echo -e ""{}"\\t0,0,178" \\; > peaks.igv.txt
    # Avoid error when consensus not produced
    find * -type l -name "*.bed" -exec echo -e ""{}"\\t0,0,178" \\; | { grep "^$consensus_dir" || test \$? = 1; } > consensus.igv.txt

    touch replace_paths.txt
    if [ -d "mappings" ]; then
        cat mappings/* > replace_paths.txt
    fi

    cat *.igv.txt > igv_files_orig.txt
    igv_files_to_session.py \\
        igv_session.xml \\
        igv_files_orig.txt \\
        replace_paths.txt \\
        ../../genome/${fasta.getName()} \\
        --path_prefix '../../' \\
        --json_out igv_session.json \\
        --genome_id '${genome_id}' \\
        ${ip_bam_args} \\
        ${ip_bai_args} \\
        ${ctrl_bam_args} \\
        ${ctrl_bai_args} \\
        ${sample_id_args} \\
        ${control_id_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
