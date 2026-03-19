//
// Uncompress and prepare reference genome files
//

include {
    GUNZIP as GUNZIP_FASTA ;
    GUNZIP as GUNZIP_GTF ;
    GUNZIP as GUNZIP_GFF ;
    GUNZIP as GUNZIP_GENE_BED ;
    GUNZIP as GUNZIP_BLACKLIST
} from '../../../modules/nf-core/gunzip/main'

include {
    UNTAR as UNTAR_BWA_INDEX ;
    UNTAR as UNTAR_BOWTIE2_INDEX ;
    UNTAR as UNTAR_STAR_INDEX
} from '../../../modules/nf-core/untar/main'

include { UNTARFILES                   } from '../../../modules/nf-core/untarfiles/main'
include { GFFREAD                      } from '../../../modules/nf-core/gffread/main'
include { SAMTOOLS_FAIDX               } from '../../../modules/nf-core/samtools/faidx/main'
include { BWA_INDEX                    } from '../../../modules/nf-core/bwa/index/main'
include { BOWTIE2_BUILD                } from '../../../modules/nf-core/bowtie2/build/main'
include { CHROMAP_INDEX                } from '../../../modules/nf-core/chromap/index/main'

include { GTF2BED                      } from '../../../modules/local/gtf2bed/main'
include { GENOME_BLACKLIST_REGIONS     } from '../../../modules/local/genome_blacklist_regions/main'
include { STAR_GENOMEGENERATE          } from '../../../modules/local/star_genomegenerate/main'

workflow PREPARE_GENOME {
    take:
    genome             //  string: genome name
    genomes            //     map: genome attributes
    prepare_tool_index //  string: tool to prepare index for
    fasta              //    path: path to genome fasta file
    gtf                //    file: /path/to/genome.gtf
    gff                //    file: /path/to/genome.gff
    blacklist          //    file: /path/to/blacklist.bed
    gene_bed           //    file: /path/to/gene.bed
    bwa_index          //    file: /path/to/bwa/index/
    bowtie2_index      //    file: /path/to/bowtie2/index/
    chromap_index      //    file: /path/to/chromap/index/
    star_index         //    file: /path/to/star/index/

    main:

    //
    // Uncompress genome fasta file if required
    //
    ch_fasta = Channel.empty()
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA([[:], fasta]).gunzip.map { it[1] }
    }
    else {
        ch_fasta = Channel.value(file(fasta, checkIfExists: true))
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (gtf) {
        if (gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF([[:], gtf]).gunzip.map { it[1] }
        } else {
            ch_gtf = Channel.value(file(gtf, checkIfExists: true))
        }
    } else if (gff) {
        if (gff.endsWith('.gz')) {
            ch_gff      = GUNZIP_GFF([[:], file(gff, checkIfExists: true)]).gunzip.map { it[1] }
        } else {
            ch_gff = Channel.value(file(gff, checkIfExists: true)).map { [ [:], it ] }
        }

        ch_gtf      = GFFREAD(ch_gff, []).gtf.map { it[1] }
    }

    //
    // Uncompress blacklist file if required
    //
    ch_blacklist = Channel.empty()
    if (blacklist) {
        if (blacklist.endsWith('.gz')) {
            ch_blacklist = GUNZIP_BLACKLIST([[:], blacklist]).gunzip.map { it[1] }
        } else {
            ch_blacklist = Channel.value(file(blacklist))
        }
    }

    //
    // Uncompress gene BED annotation file or create from GTF if required
    //

    // If --gtf is supplied along with --genome
    // Make gene bed from supplied --gtf instead of using iGenomes one automatically
    def make_bed = false
    if (!gene_bed) {
        make_bed = true
    } else if (genome && gtf) {
        if (genomes[genome].gtf != gtf) {
            make_bed = true
        }
    }

    if (make_bed) {
        ch_gene_bed = GTF2BED(ch_gtf).bed
    } else {
        if (gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED([[:], gene_bed]).gunzip.map { it[1] }
        } else {
            ch_gene_bed = Channel.value(file(gene_bed))
        }
    }

    //
    // Create chromosome sizes file
    //
    ch_chrom_sizes = channel.empty()
    ch_fai         = channel.empty()

    SAMTOOLS_FAIDX(ch_fasta.map { item -> [ [:], item, [] ] }, true)
    ch_chrom_sizes = SAMTOOLS_FAIDX.out.sizes.map { tuple -> tuple[1] }
    ch_fai         = SAMTOOLS_FAIDX.out.fai.map { tuple -> tuple[1] }

    //
    // Prepare genome intervals for filtering by removing regions in blacklist file
    //
    ch_genome_filtered_bed = Channel.empty()

    GENOME_BLACKLIST_REGIONS(
        ch_chrom_sizes,
        ch_blacklist.ifEmpty([])
    )
    ch_genome_filtered_bed = GENOME_BLACKLIST_REGIONS.out.bed

    //
    // Uncompress BWA index or generate from scratch if required
    //
    ch_bwa_index = Channel.empty()
    if (prepare_tool_index == 'bwa') {
        if (bwa_index) {
            if (bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX([[:], bwa_index]).untar
            } else {
                ch_bwa_index = [[:], file(bwa_index)]
            }
        } else {
            ch_bwa_index = BWA_INDEX(ch_fasta.map { [[:], it] }).index
        }
    }

    //
    // Uncompress Bowtie2 index or generate from scratch if required
    //
    ch_bowtie2_index = Channel.empty()
    if (prepare_tool_index == 'bowtie2') {
        if (bowtie2_index) {
            if (bowtie2_index.endsWith('.tar.gz')) {
                ch_bowtie2_index = UNTAR_BOWTIE2_INDEX([[:], bowtie2_index]).untar
            } else {
                ch_bowtie2_index = [[:], file(bowtie2_index)]
            }
        } else {
            ch_bowtie2_index = BOWTIE2_BUILD(ch_fasta.map { [[:], it] }).index
        }
    }

    //
    // Uncompress CHROMAP index or generate from scratch if required
    //
    ch_chromap_index = Channel.empty()
    if (prepare_tool_index == 'chromap') {
        if (chromap_index) {
            if (chromap_index.endsWith('.tar.gz')) {
                ch_chromap_index = UNTARFILES([[:], chromap_index]).files
            } else {
                ch_chromap_index = [[:], file(chromap_index)]
            }
        } else {
            ch_chromap_index = CHROMAP_INDEX(ch_fasta.map { [[:], it] }).index
        }
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if (prepare_tool_index == 'star') {
        if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX([[:], star_index]).untar.map { it[1] }
            } else {
                ch_star_index = Channel.value(file(star_index))
            }
        } else {
            ch_star_index = STAR_GENOMEGENERATE(ch_fasta, ch_gtf).index
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
}
