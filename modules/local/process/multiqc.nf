def SOFTWARE = 'multiqc'

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

process MULTIQC {
    label 'process_medium'
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    //container "https://depot.galaxyproject.org/singularity/multiqc:1.9--pyh9f0ad1d_0"

    conda (params.conda ? "bioconda::multiqc=1.9" : null)

    input:
    path multiqc_config
    path mqc_custom_config
    path software_versions
    val workflow_summary

    path ('fastqc/*')
    path ('trimgalore/*')
    path ('trimgalore/fastqc/*')

    // path ('alignment/library/*')
    // path ('alignment/library/*')
    // path ('alignment/library/*')

    // path ('alignment/mergedLibrary/*') from ch_merge_bam_stats_mqc.collect()
    // path ('alignment/mergedLibrary/*') from ch_rm_orphan_flagstat_mqc.collect{it[1]}
    // path ('alignment/mergedLibrary/*') from ch_rm_orphan_stats_mqc.collect()
    // path ('alignment/mergedLibrary/picard_metrics/*') from ch_merge_bam_metrics_mqc.collect()
    // path ('alignment/mergedLibrary/picard_metrics/*') from ch_collectmetrics_mqc.collect()
    //
    // path ('macs/*') from ch_macs_mqc.collect().ifEmpty([])
    // path ('macs/*') from ch_macs_qc_mqc.collect().ifEmpty([])
    // path ('macs/consensus/*') from ch_macs_consensus_counts_mqc.collect().ifEmpty([])
    // path ('macs/consensus/*') from ch_macs_consensus_deseq_mqc.collect().ifEmpty([])
    //
    // path ('preseq/*') from ch_preseq_mqc.collect().ifEmpty([])
    // path ('deeptools/*') from ch_plotfingerprint_mqc.collect().ifEmpty([])
    // path ('deeptools/*') from ch_plotprofile_mqc.collect().ifEmpty([])
    // path ('phantompeakqualtools/*') from ch_spp_out_mqc.collect().ifEmpty([])
    // path ('phantompeakqualtools/*') from ch_spp_csv_mqc.collect().ifEmpty([])

    val opts

    output:
    path "*multiqc_report.html", emit: report
    path "*_data", emit: data

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    // TODO nf-core: Specify which MultiQC modules to use with -m for a faster run time
    """
    echo '$workflow_summary' > workflow_summary_mqc.yaml
    multiqc -f $opts.args $rtitle $rfilename $custom_config_file .
    multiqc --version | sed -e "s/multiqc, version //g" > ${SOFTWARE}.version.txt
    """
}

// path ('fastqc/*') from ch_fastqc_reports_mqc.collect().ifEmpty([])
// path ('trimgalore/*') from ch_trimgalore_results_mqc.collect().ifEmpty([])
// path ('trimgalore/fastqc/*') from ch_trimgalore_fastqc_reports_mqc.collect().ifEmpty([])
//
// path ('alignment/library/*') from ch_sort_bam_flagstat_mqc.collect()
// path ('alignment/mergedLibrary/*') from ch_merge_bam_stats_mqc.collect()
// path ('alignment/mergedLibrary/*') from ch_rm_orphan_flagstat_mqc.collect{it[1]}
// path ('alignment/mergedLibrary/*') from ch_rm_orphan_stats_mqc.collect()
// path ('alignment/mergedLibrary/picard_metrics/*') from ch_merge_bam_metrics_mqc.collect()
// path ('alignment/mergedLibrary/picard_metrics/*') from ch_collectmetrics_mqc.collect()
//
// path ('macs/*') from ch_macs_mqc.collect().ifEmpty([])
// path ('macs/*') from ch_macs_qc_mqc.collect().ifEmpty([])
// path ('macs/consensus/*') from ch_macs_consensus_counts_mqc.collect().ifEmpty([])
// path ('macs/consensus/*') from ch_macs_consensus_deseq_mqc.collect().ifEmpty([])
//
// path ('preseq/*') from ch_preseq_mqc.collect().ifEmpty([])
// path ('deeptools/*') from ch_plotfingerprint_mqc.collect().ifEmpty([])
// path ('deeptools/*') from ch_plotprofile_mqc.collect().ifEmpty([])
// path ('phantompeakqualtools/*') from ch_spp_out_mqc.collect().ifEmpty([])
// path ('phantompeakqualtools/*') from ch_spp_csv_mqc.collect().ifEmpty([])
