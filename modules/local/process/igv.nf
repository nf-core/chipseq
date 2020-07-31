/*
 * Create IGV session file
 */
process IGV {
    publishDir "${params.outdir}/${options.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (options.publish_results == "none") null
                      else filename }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path fasta
    path ("${bigwig_options.publish_dir}/*")
    path ("${peak_options.publish_dir}/*")
    path ("${consensus_options.publish_dir}/*")
    val bigwig_options
    val peak_options
    val consensus_options
    // path differential_peaks from ch_macs_consensus_deseq_comp_igv.collect().ifEmpty([])
    val options

    output:
    path "*files.txt", emit: txt
    path "*.xml", emit: xml

    script: // scripts are bundled with the pipeline in nf-core/chipseq/bin/
    """
    find * -type l -name "*.bigWig" -exec echo -e ""{}"\\t0,0,178" \\; > bigwig.igv.txt
    find * -type l -name "*Peak" -exec echo -e ""{}"\\t0,0,178" \\; > peaks.igv.txt

    cat *.txt > igv_files.txt
    igv_files_to_session.py igv_session.xml igv_files.txt ../../genome/${fasta.getName()} --path_prefix '../../'
    """
}
// find * -type f -name "${prefix}.bed" -exec echo -e "bwa/mergedLibrary/macs/${PEAK_TYPE}/consensus/${antibody}/"{}"\\t0,0,0" \\; > ${prefix}.bed.igv.txt
// find * -type f -name "*.FDR0.05.results.bed" -exec echo -e "bwa/mergedLibrary/macs/${PEAK_TYPE}/consensus/${antibody}/deseq2/"{}"\\t255,0,0" \\; > ${prefix}.igv.txt
