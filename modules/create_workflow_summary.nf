// Initialise parameters
summary = null

/*
 * Create summary file for MultiQC
 */
def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-chipseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/chipseq Workflow Summary'
    section_href: 'https://github.com/nf-core/chipseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}
