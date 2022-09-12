#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/chipseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/chipseq
    Website: https://nf-co.re/chipseq
    Slack  : https://nfcore.slack.com/channels/chipseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta         = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.bwa_index     = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.bowtie2_index = WorkflowMain.getGenomeAttribute(params, 'bowtie2')
params.chromap_index = WorkflowMain.getGenomeAttribute(params, 'chromap')
params.star_index    = WorkflowMain.getGenomeAttribute(params, 'star')
params.gtf           = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff           = WorkflowMain.getGenomeAttribute(params, 'gff')
params.gene_bed      = WorkflowMain.getGenomeAttribute(params, 'gene_bed')
params.blacklist     = WorkflowMain.getGenomeAttribute(params, 'blacklist')
params.macs_gsize    = WorkflowMain.getMacsGsize(params)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHIPSEQ } from './workflows/chipseq'

//
// WORKFLOW: Run main nf-core/chipseq analysis pipeline
//
workflow NFCORE_CHIPSEQ {
    CHIPSEQ ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_CHIPSEQ ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
