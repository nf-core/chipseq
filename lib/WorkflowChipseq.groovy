//
// This file holds several functions specific to the workflow/chipseq.nf in the nf-core/chipseq pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowChipseq {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, valid_params) {
        genomeExistsError(params, log)


        if (!params.fasta) {
            Nextflow.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
        }

        if (!params.gtf && !params.gff) {
            def error_string = "No GTF or GFF3 annotation specified! The pipeline requires at least one of these files."
            Nextflow.error(error_string)
        }

        if (params.gtf && params.gff) {
            gtfGffWarn(log)
        }

        if (!params.macs_gsize) {
            macsGsizeWarn(log)
        }

        if (!params.read_length && !params.macs_gsize) {
            def error_string = "Both '--read_length' and '--macs_gsize' not specified! Please specify either to infer MACS2 genome size for peak calling."
            Nextflow.error(error_string)
        }

        if (params.aligner) {
            if (!valid_params['aligners'].contains(params.aligner)) {
                    def error_string = "Invalid option: '${params.aligner}'. Valid options for '--aligner': ${valid_params['aligners'].join(', ')}."
                    Nextflow.error(error_string)
            }
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }

    //
    // Print a warning if both GTF and GFF have been provided
    //
    private static void gtfGffWarn(log) {
        log.warn "=============================================================================\n" +
            "  Both '--gtf' and '--gff' parameters have been provided.\n" +
            "  Using GTF file as priority.\n" +
            "==================================================================================="
    }

    //
    // Print a warning if macs_gsize parameter has not been provided
    //
    private static void macsGsizeWarn(log) {
        log.warn "=============================================================================\n" +
            "  --macs_gsize parameter has not been provided.\n" +
            "  It will be auto-calculated by 'khmer unique-kmers.py' using the '--read_length' parameter.\n" +
            "  Explicitly provide '--macs_gsize' to change this behaviour.\n" +
            "==================================================================================="
    }

}
