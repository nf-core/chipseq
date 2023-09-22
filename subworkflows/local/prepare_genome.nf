//
// Uncompress and prepare reference genome files
//

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_GENE_BED
    GUNZIP as GUNZIP_BLACKLIST } from '../../modules/nf-core/gunzip/main'

include {
    UNTAR as UNTAR_BWA_INDEX
    UNTAR as UNTAR_BOWTIE2_INDEX
    UNTAR as UNTAR_STAR_INDEX    } from '../../modules/nf-core/untar/main'

include { UNTARFILES           } from '../../modules/nf-core/untarfiles/main'
include { GFFREAD              } from '../../modules/nf-core/gffread/main'
include { CUSTOM_GETCHROMSIZES } from '../../modules/nf-core/custom/getchromsizes/main'
include { BWA_INDEX            } from '../../modules/nf-core/bwa/index/main'
include { BOWTIE2_BUILD        } from '../../modules/nf-core/bowtie2/build/main'
include { CHROMAP_INDEX        } from '../../modules/nf-core/chromap/index/main'

include { GTF2BED                  } from '../../modules/local/gtf2bed'
include { GENOME_BLACKLIST_REGIONS } from '../../modules/local/genome_blacklist_regions'
include { STAR_GENOMEGENERATE      } from '../../modules/local/star_genomegenerate'

workflow PREPARE_GENOME {
    take:
    prepare_tool_index // string  : tool to prepare index for

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    ch_fasta = Channel.empty()
    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.map{ it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(params.fasta))
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], params.gtf ] ).gunzip.map{ it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = file(params.gtf)
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff      = GUNZIP_GFF ( [ [:], params.gff ] ).gunzip.map{ it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = file(params.gff)
        }
        ch_gtf      = GFFREAD ( ch_gff ).gtf
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

    //
    // Uncompress blacklist file if required
    //
    ch_blacklist = Channel.empty()
    if (params.blacklist) {
        if (params.blacklist.endsWith('.gz')) {
            ch_blacklist = GUNZIP_BLACKLIST ( [ [:], params.blacklist ] ).gunzip.map{ it[1] }
            ch_versions  = ch_versions.mix(GUNZIP_BLACKLIST.out.versions)
        } else {
            ch_blacklist = Channel.fromPath(file(params.blacklist))
        }
    }

    //
    // Uncompress gene BED annotation file or create from GTF if required
    //

    // If --gtf is supplied along with --genome
    // Make gene bed from supplied --gtf instead of using iGenomes one automatically
    def make_bed = false
    if (!params.gene_bed) {
        make_bed = true
    } else if (params.genome && params.gtf) {
        if (params.genomes[ params.genome ].gtf != params.gtf) {
            make_bed = true
        }
    }

    if (make_bed) {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    } else {
        if (params.gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( [ [:], params.gene_bed ] ).gunzip.map{ it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = Channel.value(file(params.gene_bed))
        }
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map{ it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Prepare genome intervals for filtering by removing regions in blacklist file
    //
    ch_genome_filtered_bed = Channel.empty()

    GENOME_BLACKLIST_REGIONS (
        CUSTOM_GETCHROMSIZES.out.sizes.map{ it[1] },
        ch_blacklist.ifEmpty([])
    )
    ch_genome_filtered_bed = GENOME_BLACKLIST_REGIONS.out.bed
    ch_versions = ch_versions.mix(GENOME_BLACKLIST_REGIONS.out.versions)


    //
    // Uncompress BWA index or generate from scratch if required
    //
    ch_bwa_index = Channel.empty()
    if (prepare_tool_index == 'bwa') {
        if (params.bwa_index) {
            if (params.bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( [ [:], params.bwa_index ] ).untar
                ch_versions  = ch_versions.mix(UNTAR_BWA_INDEX.out.versions)
            } else {
                ch_bwa_index = file(params.bwa_index)
            }
        } else {
            ch_bwa_index = BWA_INDEX ( ch_fasta.map { [ [:], it ] } ).index
            ch_versions  = ch_versions.mix(BWA_INDEX.out.versions)
        }
    }

    //
    // Uncompress Bowtie2 index or generate from scratch if required
    //
    ch_bowtie2_index = Channel.empty()
    if (prepare_tool_index == 'bowtie2') {
        if (params.bowtie2_index) {
            if (params.bowtie2_index.endsWith('.tar.gz')) {
                ch_bowtie2_index = UNTAR_BOWTIE2_INDEX ( [ [:], params.bowtie2_index ] ).untar
                ch_versions  = ch_versions.mix(UNTAR_BOWTIE2_INDEX.out.versions)
            } else {
                ch_bowtie2_index = [ [:], file(params.bowtie2_index) ]
            }
        } else {
            ch_bowtie2_index = BOWTIE2_BUILD ( ch_fasta.map { [ [:], it ] } ).index
            ch_versions      = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }
    }

    //
    // Uncompress CHROMAP index or generate from scratch if required
    //
    ch_chromap_index = Channel.empty()
    if (prepare_tool_index == 'chromap') {
        if (params.chromap_index) {
            if (params.chromap_index.endsWith('.tar.gz')) {
                ch_chromap_index = UNTARFILES ( [ [:], params.chromap_index ] ).files
                ch_versions  = ch_versions.mix(UNTARFILES.out.versions)
            } else {
                ch_chromap_index = [ [:], file(params.chromap_index) ]
            }
        } else {
            ch_chromap_index = CHROMAP_INDEX ( ch_fasta.map { [ [:], it ] } ).index
            ch_versions  = ch_versions.mix(CHROMAP_INDEX.out.versions)
        }
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if (prepare_tool_index == 'star') {
        if (params.star_index) {
            if (params.star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( [ [:], params.star_index ] ).untar.map{ it[1] }
                ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = Channel.value(file(params.star_index))
            }
        } else {
            ch_star_index = STAR_GENOMEGENERATE ( ch_fasta, ch_gtf ).index
            ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        }
    }

    emit:
    fasta         = ch_fasta                  //    path: genome.fasta
    fai           = ch_fai                    //    path: genome.fai
    gtf           = ch_gtf                    //    path: genome.gtf
    gene_bed      = ch_gene_bed               //    path: gene.bed
    chrom_sizes   = ch_chrom_sizes            //    path: genome.sizes
    filtered_bed  = ch_genome_filtered_bed    //    path: *.include_regions.bed
    bwa_index     = ch_bwa_index              //    path: bwa/index/
    bowtie2_index = ch_bowtie2_index          //    path: bowtie2/index/
    chromap_index = ch_chromap_index          //    path: genome.index
    star_index    = ch_star_index             //    path: star/index/

    versions    = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
