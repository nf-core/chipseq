// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FRIP_SCORE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.conda ? "${baseDir}/environment.yml" : null) //TODO update with mulled env (bedtools+samtools)

    input:
    tuple val(meta), path(bam), path(peak)

    output:
    tuple val(meta), path("*.txt"), emit: txt

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    READS_IN_PEAKS=\$(intersectBed -a $bam -b $peak $options.args | awk -F '\t' '{sum += \$NF} END {print sum}')
    samtools flagstat $bam > ${bam}.flagstat
    grep 'mapped (' ${bam}.flagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${prefix}", a/\$1}' > ${prefix}.FRiP.txt
    """
}
