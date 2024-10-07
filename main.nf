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

params.fasta         = getGenomeAttribute('fasta')
params.bwa_index     = getGenomeAttribute('bwa')
params.bowtie2_index = getGenomeAttribute('bowtie2')
params.chromap_index = getGenomeAttribute('chromap')
params.star_index    = getGenomeAttribute('star')
params.gtf           = getGenomeAttribute('gtf')
params.gff           = getGenomeAttribute('gff')
params.gene_bed      = getGenomeAttribute('gene_bed')
params.blacklist     = getGenomeAttribute('blacklist')
params.macs_gsize    = getMacsGsize(params)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CHIPSEQ                 } from './workflows/chipseq'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_chipseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_chipseq_pipeline'
// include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_chipseq_pipeline'
// include { getMacsGsize            } from './subworkflows/local/utils_nfcore_chipseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_CHIPSEQ {

    main:
    ch_versions = Channel.empty()

    // SUBWORKFLOW: Prepare genome files
    PREPARE_GENOME (
        params.genome,
        params.genomes,
        params.aligner,
        params.fasta,
        params.gtf,
        params.gff,
        params.blacklist,
        params.gene_bed,
        params.bwa_index,
        params.bowtie2_index,
        params.chromap_index,
        params.star_index,
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // WORKFLOW: Run nf-core/chipseq workflow
    //
    ch_input = Channel.value(file(params.input, checkIfExists: true))

    CHIPSEQ(
        ch_input,
        ch_versions,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.gene_bed,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.filtered_bed,
        PREPARE_GENOME.out.bwa_index,
        PREPARE_GENOME.out.bowtie2_index,
        PREPARE_GENOME.out.chromap_index,
        PREPARE_GENOME.out.star_index
    )

    emit:
    multiqc_report = CHIPSEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions                // channel: [version1, version2, ...]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_CHIPSEQ ()

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_CHIPSEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Get macs genome size (macs_gsize)
//
def getMacsGsize(params) {
    def val = null
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey('macs_gsize')) {
            if (params.genomes[ params.genome ][ 'macs_gsize' ].containsKey(params.read_length.toString())) {
                val = params.genomes[ params.genome ][ 'macs_gsize' ][ params.read_length.toString() ]
            }
        }
    }
    return val
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
