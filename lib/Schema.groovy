/*
 * This file holds several functions used to perform JSON parameter validation, help and summary rendering for the nf-core pipeline template.
 */

import groovy.json.JsonSlurper

class Schema {
    /*
     * This method tries to read a JSON params file
     */
    private static LinkedHashMap params_get(String path) {
        def params_map = new LinkedHashMap()
        try {
            params_map = params_try(path)
        } catch (Exception e) {
            println "Could not read parameters settings from JSON. $e"
            params_map = new LinkedHashMap()
        }
        return params_map
    }

    /*
    Method to actually read in JSON file using Groovy.
    Group (as Key), values are all parameters
        - Parameter1 as Key, Description as Value
        - Parameter2 as Key, Description as Value
        ....
    Group
        -
    */
    private static LinkedHashMap params_try(String path) throws Exception {

        def json = new File(path).text
        def Map json_params = (Map) new JsonSlurper().parseText(json).get('properties')

        /* Tree looks like this in nf-core schema
        *  properties <- this is what the first get('properties') gets us
            group 1
                properties
                description
            group 2
                properties
                description
            group 3
                properties
                description
        */
        def params_map = new LinkedHashMap()
        json_params.each { key, val ->
            def Map group = json_params."$key".properties // Gets the property object of the group
            def sub_params = new LinkedHashMap()
            group.each { innerkey, value ->
                sub_params.put("$innerkey", [ "$value.type", "$value.description" ])
            }
            params_map.put("$key", sub_params)
        }
        return params_map
    }

    private static Integer params_max_chars(params_map) {
        Integer max_chars = 0
        for (group in params_map.keySet()) {
            def params = params_map.get(group)  // This gets the parameters of that particular group
            for (par in params.keySet()) {
                if (par.size() > max_chars) {
                    max_chars = par.size()
                }
            }
        }
        return max_chars
    }

    private static String params_beautify(params_map) {
        String output = ""
        def max_chars = params_max_chars(params_map) + 1
        for (group in params_map.keySet()) {
            output += group + "\n"
            def params = params_map.get(group)  // This gets the parameters of that particular group
            for (par in params.keySet()) {
                def type = "[" + params.get(par)[0] + "]"
                def description = params.get(par)[1]
                output+= "    \u001B[1m" +  par.padRight(max_chars) + "\u001B[1m" + type.padRight(10) + description + "\n"
            }
            output += "\n"
        }
        return output
    }

    private static String params_help(path, command) {
        String output = "Typical pipeline command:\n\n"
        output += "    ${command}\n\n"
        output += params_beautify(params_get(path))
    }

    private static LinkedHashMap params_summary(workflow, params, run_name) {
        def Map summary = [:]
        if (workflow.revision) summary['Pipeline Release'] = workflow.revision
        summary['Run Name']         = run_name ?: workflow.runName
        summary['Design File']            = params.input
        summary['Genome']                 = params.genome ?: 'Not supplied'
        summary['Fasta File']             = params.fasta
        summary['GTF File']               = params.gtf
        if (params.gene_bed)              summary['Gene BED File'] = params.gene_bed
        if (params.bwa_index)             summary['BWA Index'] = params.bwa_index
        if (params.blacklist)             summary['Blacklist BED'] = params.blacklist
        if (params.bwa_min_score)         summary['BWA Min Score'] = params.bwa_min_score
        summary['MACS2 Genome Size']      = params.macs_gsize ?: 'Not supplied'
        summary['Min Consensus Reps']     = params.min_reps_consensus
        if (params.macs_gsize)            summary['MACS2 Narrow Peaks'] = params.narrow_peak ? 'Yes' : 'No'
        if (!params.narrow_peak)          summary['MACS2 Broad Cutoff'] = params.broad_cutoff
        if (params.macs_fdr)              summary['MACS2 FDR'] = params.macs_fdr
        if (params.macs_pvalue)           summary['MACS2 P-value'] = params.macs_pvalue
        if (params.skip_trimming) {
            summary['Trimming Step']      = 'Skipped'
        } else {
            summary['Trim R1']            = "$params.clip_r1 bp"
            summary['Trim R2']            = "$params.clip_r2 bp"
            summary["Trim 3' R1"]         = "$params.three_prime_clip_r1 bp"
            summary["Trim 3' R2"]         = "$params.three_prime_clip_r2 bp"
            summary['NextSeq Trim']       = "$params.trim_nextseq bp"
        }
        if (params.seq_center)            summary['Sequencing Center'] = params.seq_center
        summary['Fragment Size']          = "$params.fragment_size bp"
        summary['Fingerprint Bins']       = params.fingerprint_bins
        if (params.keep_dups)             summary['Keep Duplicates'] = 'Yes'
        if (params.keep_multi_map)        summary['Keep Multi-mapped'] = 'Yes'
        summary['Save Genome Index']      = params.save_reference ? 'Yes' : 'No'
        if (params.save_trimmed)          summary['Save Trimmed'] = 'Yes'
        if (params.save_align_intermeds)  summary['Save Intermeds'] =  'Yes'
        if (params.save_macs_pileup)      summary['Save MACS2 Pileup'] = 'Yes'
        if (params.skip_peak_qc)          summary['Skip MACS2 Peak QC'] = 'Yes'
        if (params.skip_peak_annotation)  summary['Skip Peak Annotation'] = 'Yes'
        if (params.skip_consensus_peaks)  summary['Skip Consensus Peaks'] = 'Yes'
        if (params.deseq2_vst)            summary['Use DESeq2 vst Transform'] = 'Yes'
        if (params.skip_diff_analysis)    summary['Skip Differential Analysis'] = 'Yes'
        if (params.skip_fastqc)           summary['Skip FastQC'] = 'Yes'
        if (params.skip_picard_metrics)   summary['Skip Picard Metrics'] = 'Yes'
        if (params.skip_preseq)           summary['Skip Preseq'] = 'Yes'
        if (params.skip_plot_profile)     summary['Skip plotProfile'] = 'Yes'
        if (params.skip_plot_fingerprint) summary['Skip plotFingerprint'] = 'Yes'
        if (params.skip_spp)              summary['Skip spp'] = 'Yes'
        if (params.skip_igv)              summary['Skip IGV'] = 'Yes'
        if (params.skip_multiqc)          summary['Skip MultiQC'] = 'Yes'
        summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
        if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
        summary['Output dir']       = params.outdir
        summary['Launch dir']       = workflow.launchDir
        summary['Working dir']      = workflow.workDir
        summary['Script dir']       = workflow.projectDir
        summary['User']             = workflow.userName
        if (workflow.profile.contains('awsbatch')) {
            summary['AWS Region']   = params.awsregion
            summary['AWS Queue']    = params.awsqueue
            summary['AWS CLI']      = params.awscli
        }
        summary['Config Profile'] = workflow.profile
        if (params.config_profile_description) summary['Config Profile Descr']   = params.config_profile_description
        if (params.config_profile_contact)     summary['Config Profile Contact'] = params.config_profile_contact
        if (params.config_profile_url)         summary['Config Profile URL']     = params.config_profile_url
        summary['Config Files'] = workflow.configFiles.join(', ')
        if (params.email || params.email_on_fail) {
            summary['E-mail Address']    = params.email
            summary['E-mail on failure'] = params.email_on_fail
            summary['MultiQC maxsize']   = params.max_multiqc_email_size
        }
        return summary
    }

    static String params_mqc_summary(summary) {
        String yaml_file_text  = """
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

        return yaml_file_text
    }
}
