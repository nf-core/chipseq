def SOFTWARE = 'macs2'

process MACS2_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else if (filename.endsWith('.version.txt')) null
                      else filename }

    container "quay.io/biocontainers/macs2:2.2.7.1--py37h516909a_0"
    //container "https://depot.galaxyproject.org/singularity/macs2:2.2.7.1--py37h516909a_0"

    conda (params.conda ? "bioconda::macs2=2.2.7.1" : null)

    input:
    tuple val(meta), path(ipbam), path(controlbam)
    val macs2_gsize
    val opts

    output:
    tuple val(meta), path("*.{narrowPeak,broadPeak}"), emit: peak
    path '*.xls', emit: xls
    path '*.gappedPeak', emit: gapped optional true
    path '*.bed', emit: bed optional true
    path '*.bdg', emit: bdg optional true
    path "*.version.txt", emit: version

    script:
    prefix = opts.suffix ? "${meta.id}${opts.suffix}" : "${meta.id}"
    format = meta.single_end ? 'BAM' : 'BAMPE'
    control = controlbam ? "--control $controlbam" : ''
    """
    macs2 \\
        callpeak \\
        $opts.args \\
        --gsize $macs2_gsize \\
        --format $format \\
        --name $prefix \\
        --treatment $ipbam \\
         $control

    macs2 --version | sed -e "s/macs2 //g" > ${SOFTWARE}.version.txt
    """
}
// cat ${ip}_peaks.${PEAK_TYPE} | wc -l | awk -v OFS='\t' '{ print "${ip}", \$1 }' | cat $peak_count_header - > ${ip}_peaks.count_mqc.tsv
//
// find * -type f -name "*.${PEAK_TYPE}" -exec echo -e "bwa/mergedLibrary/macs/${PEAK_TYPE}/"{}"\\t0,0,178" \\; > ${ip}_peaks.igv.txt
