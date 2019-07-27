// Initialise parameters
params.outdir = './results'

process check_design {
    tag "$design"
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file design

    output:
    file "design_reads.csv"
    file "design_controls.csv"

    script:  // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    """
    check_design.py $design design_reads.csv design_controls.csv
    """
}
