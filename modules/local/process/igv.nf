// /*
//  * STEP 8: Create IGV session file
//  */
// process IGV {
//     publishDir "${params.outdir}/igv/${PEAK_TYPE}", mode: params.publish_dir_mode
//
//     when:
//     !params.skip_igv
//
//     input:
//     path fasta from ch_fasta
//     path bigwigs from ch_bigwig_igv.collect().ifEmpty([])
//     path peaks from ch_macs_igv.collect().ifEmpty([])
//     path consensus_peaks from ch_macs_consensus_igv.collect().ifEmpty([])
//     path differential_peaks from ch_macs_consensus_deseq_comp_igv.collect().ifEmpty([])
//
//     output:
//     path '*.{txt,xml}'
//
//     script: // scripts are bundled with the pipeline in nf-core/chipseq/bin/
//     """
//     cat *.txt > igv_files.txt
//     igv_files_to_session.py igv_session.xml igv_files.txt ../../genome/${fasta.getName()} --path_prefix '../../'
//     """
// }
