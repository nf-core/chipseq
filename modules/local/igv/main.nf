/*
 * Create IGV session file (XML for IGV Desktop, JSON for IGV.js / Seqera Data Explorer)
 */
process IGV {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/62/622d8944750bc95bb56b4c3ed5c2b827e677c14073d48a5231e0f2bec0718add/data' :
        'community.wave.seqera.io/library/python:3.12.12--74abbf3898230efd' }"

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
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script: // scripts are bundled with the pipeline in nf-core/chipseq/bin/
    def consensus_dir = "${aligner_dir}/merged_library/macs3/${peak_dir}/consensus/*"

    // Build space-separated file lists for the Python script.
    // Inputs may be single files, lists of files, or empty lists [].
    def ip_bam_list     = ip_bams instanceof List ? ip_bams : (ip_bams.name != 'input' ? [ip_bams] : [])
    def ip_bai_list     = ip_bais instanceof List ? ip_bais : (ip_bais.name != 'input' ? [ip_bais] : [])
    def ctrl_bam_list   = control_bams instanceof List ? control_bams : (control_bams.name != 'input' ? [control_bams] : [])
    def ctrl_bai_list   = control_bais instanceof List ? control_bais : (control_bais.name != 'input' ? [control_bais] : [])
    def ip_bam_args     = ip_bam_list   ? "--bam_files ${ip_bam_list.collect{ it.name }.join(' ')}"           : ''
    def ip_bai_args     = ip_bai_list   ? "--bai_files ${ip_bai_list.collect{ it.name }.join(' ')}"           : ''
    def ctrl_bam_args   = ctrl_bam_list ? "--control_bam_files ${ctrl_bam_list.collect{ it.name }.join(' ')}" : ''
    def ctrl_bai_args   = ctrl_bai_list ? "--control_bai_files ${ctrl_bai_list.collect{ it.name }.join(' ')}" : ''
    def sample_id_args  = ip_sample_ids     ? "--sample_ids ${ip_sample_ids.join(' ')}"      : ''
    def control_id_args = control_sample_ids ? "--control_ids ${control_sample_ids.join(' ')}" : ''
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
    """
}
