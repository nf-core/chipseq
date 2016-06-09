#!/usr/bin/env nextflow

/*
========================================================================================
                    ChIP-seq    v1.0
========================================================================================
 ChIP-seq Analysis Pipeline.
 @Authors
 Chuan Wang <chuan.wang@scilifelab.se>
 Release 2016-06-08
----------------------------------------------------------------------------------------
 Basic command:
 $ nextflow chipseq.nf

 Pipeline variables can be configured with the following command line options:
 --genome [GRCh37 | GRCm38 | NCBIM37]
 --reads [path to input files]
 --mode [single | paired]
 --macsconfig [path to config file for MACS: line format: ChIPSampleID,CtrlSampleID,AnalysisID]

 For example:
 $ nextflow main.nf -c ~/.nextflow/config --reads '*.fastq' --genome GRCh37 --mode single --macsconfig 'macssetup.config'
 $ nextflow main.nf -c ~/.nextflow/config --reads '*.R{1,2}.fastq' --genome GRCh37 --mode paired --macsconfig 'macssetup.config'
---------------------------------------------------------------------------------------
 The pipeline can determine whether the input data is single or paired end. This relies on
 specifying the input files correctly. For paired en data us the example above, i.e.
 'sample_*_{1,2}.fastq.gz'. Without the glob {1,2} (or similiar) the data will be treated
 as single end.
----------------------------------------------------------------------------------------
 Pipeline overview:
 - FastQC
 - TrimGalore
 - bwa
 - picard
 - phantompeakqualtools
 - deeptools
 - ngsplot
 - macs2
 - MultiQC
----------------------------------------------------------------------------------------
*/



/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 1.0

// Reference genome index
params.genome = 'GRCh37'
params.genomebwa = params.genomes[ params.genome ].bwa

single='null'

params.name = "ChIP-seq"

// Input files
params.reads = "data/*.fastq"
params.macsconfig = "data/*.config"

if (!params.genome) {
  exit 1, "Please specify the genome:'GRCh37','NCBIM37',or 'GRCm38'"
}

if (!params.macsconfig) {
  exit 1, "Please specify the config file for MACS"
}

// Output path
params.out = "$PWD"

// R library locations
params.rlocation = "$HOME/R/nxtflow_libs/"
nxtflow_libs=file(params.rlocation)

log.info "===================================="
log.info " ChIP-seq: v${version}"
log.info "===================================="
log.info "Reads        : ${params.reads}"
log.info "Genome       : ${params.genome}"
log.info "Genomebwa    : ${params.genomebwa}"
log.info "Current home : $HOME"
log.info "Current user : $USER"
log.info "Current path : $PWD"
log.info "R libraries  : ${params.rlocation}"
log.info "Script dir   : $baseDir"
log.info "Working dir  : $workDir"
log.info "Output dir   : ${params.out}"
log.info "===================================="

// Create R library directories if not already existing
nxtflow_libs.mkdirs()

// Set up nextflow objects
genomebwa=file(params.genomebwa)
extendReadsLen=100
macsconfig=file(params.macsconfig)

// Validate inputs

//Setting up a directory to save results to
results_path = './results'

/*
 * Create a channel for read files
 */

Channel
     .fromPath( params.reads )
     .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
     .map { path ->
        def prefix = readPrefix(path, params.reads)
        tuple(prefix, path)
     }
     .groupTuple(sort: true)
     .set { read_files }

read_files.into  { read_files_fastqc; read_files_trimming; name_for_bwa; name_for_picard; name_for_phantompeakqualtools }

/*
 * Create a channel for macs config file
 */

 MACSpara = Channel
 .from(macsconfig.readLines())
 .map { line ->
   list = line.split(',')
   ChIPSampleID = list[0]
   CtrlSampleID = list[1]
   AnalysisID = list[2]
   [ ChIPSampleID, CtrlSampleID, AnalysisID ]
 }

/*
 * STEP 1 - FastQC
 */

process fastqc {
    tag "$name"

    module 'bioinfo-tools'
    module 'FastQC'

    cpus 1
    memory '4 GB'
    time '12h'

    publishDir "$results_path/fastqc"

    input:
    set val(name), file(reads:'*') from read_files_fastqc

    output:
    file '*_fastqc.html' into fastqc_html
    file '*_fastqc.zip' into results

    """
    fastqc -q ${reads}
    """
}


/*
 * STEP 2 - Trim Galore!
 */

process trim_galore {
    tag "$name"

    module 'bioinfo-tools'
    module 'FastQC'
    module 'cutadapt'
    module 'TrimGalore'

    cpus 1
    memory '4 GB'
    time '24h'

    publishDir "$results_path/trim_galore"

    input:
    set val(name), file(reads:'*') from read_files_trimming


    output:
    file '*fq' into trimmed_reads
    file '*trimming_report.txt' into trim_galore_report

    script:
    single = reads instanceof Path
    if( !single ) {

        """
        trim_galore --paired $reads
        """

    }
    else {
        """
        trim_galore $reads
        """
    }
}

 /*
 * STEP 3 - align with bwa
 */

process bwa {
    tag "$prefix"

    module 'bioinfo-tools'
    module 'bwa'
    module 'bamtools'
    module 'samtools'
    module 'BEDTools'

    cpus 2
    memory '16 GB'
    time '120h'

    publishDir "$results_path/bwa"

    input:
    file (reads:'*') from trimmed_reads
    set val(prefix) from name_for_bwa

    output:
    file '*.sort.bam' into bam_picard
    file '*.sort.bam.bai' into results
    file '*.sort.bed' into results

    """
    bwa mem -M $genomebwa ${reads} > ${prefix}.sam
    samtools view -bT $genomebwa ${prefix}.sam > ${prefix}.bam
    samtools sort ${prefix}.bam ${prefix}.sort
    samtools index ${prefix}.sort.bam
    rm ${prefix}.sam
    rm ${prefix}.bam
    bedtools bamtobed -i ${prefix}.sort.bam | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 > ${prefix}.sort.bed
    """

}


/*
* STEP 4 Picard
*/

process picard {
    tag "$bam_picard"

    module 'bioinfo-tools'
    module 'picard'
    module 'BEDTools'
    module 'bamtools'
    module 'samtools'
    module 'java/sun_jdk1.8.0_40'

    cpus 2
    memory '16 GB'
    time '24h'

    publishDir "$results_path/picard"
    input:
    file bam_picard
    set val(prefix) from name_for_picard

    output:
    file '*.dedup.sort.bam' into bam_dedup_phantompeakqualtools,bam_dedup_ngsplotconfig,bam_dedup_ngsplot,bam_dedup_deepTools,bam_dedup_macs
    file '*.dedup.sort.bam.bai' into bai_dedup_deepTools,bai_dedup_ngsplot,bai_dedup_macs
    file '*.dedup.sort.bed' into results
    file '*.picardDupMetrics.txt' into picard_report


    """
    java -Xmx2g -jar /sw/apps/bioinfo/picard/2.0.1/milou/picard.jar MarkDuplicates INPUT=${bam_picard} OUTPUT=${prefix}.dedup.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=${prefix}.picardDupMetrics.txt VALIDATION_STRINGENCY=LENIENT PROGRAM_RECORD_ID='null'

    samtools sort ${prefix}.dedup.bam ${prefix}.dedup.sort
    samtools index ${prefix}.dedup.sort.bam
    bedtools bamtobed -i ${prefix}.dedup.sort.bam | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 > ${prefix}.dedup.sort.bed
    rm ${prefix}.dedup.bam
    """
}

/*
* STEP 5 Phantompeakqualtools
*/

process phantompeakqualtools {
    tag "$bam_dedup_phantompeakqualtools"

    module 'bioinfo-tools'
    module 'samtools'
    module 'R/3.1.0'
    module 'phantompeakqualtools'
    
    cpus 2
    memory '16 GB'
    time '24h'

    publishDir "$results_path/phantompeakqualtools"
    input:
    file bam_dedup_phantompeakqualtools
    set val(prefix) from name_for_phantompeakqualtools

    output:
    file '*.pdf' into results
    file '*.spp.out' into results

    """
    run_spp.R -c=${bam_dedup_phantompeakqualtools} -savp -out=${prefix}.spp.out
    """
}


/*
* STEP 6 deepTools
*/

process deepTools {

   module 'bioinfo-tools'
   module 'deepTools'

   cpus 4
   memory '32 GB'
   time '120h'

   publishDir "$results_path/deepTools"

   input:
   file bam from bam_dedup_deepTools.toSortedList()
   file bai from bai_dedup_deepTools.toSortedList()

   output:
   file 'fingerprints.pdf' into results
   file 'multiBamSummary.npz' into results
   file 'scatterplot_PearsonCorr_multiBamSummary.png' into results
   file 'heatmap_SpearmanCorr_multiBamSummary.png' into results

   """
   plotFingerprint -b ${bam} --plotFile fingerprints.pdf --extendReads=${extendReadsLen} --skipZeros --ignoreDuplicates --numberOfSamples 50000 --binSize=500 --plotFileFormat=pdf --plotTitle="Fingerprints"

   multiBamSummary bins --binSize=10000 --extendReads=${extendReadsLen} --ignoreDuplicates --centerReads --bamfiles ${bam} -out multiBamSummary.npz

   plotCorrelation -in multiBamSummary.npz --corMethod pearson --skipZeros --removeOutliers --plotTitle "Pearson Correlation of Read Counts" --whatToPlot scatterplot -o scatterplot_PearsonCorr_multiBamSummary.png

   plotCorrelation -in multiBamSummary.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o heatmap_SpearmanCorr_multiBamSummary.png

   """
}


/*
* STEP 7 Generate config file for ngsplot
*/

process ngs_config_generate {

    cpus 1
    memory '2 GB'
    time '1h'

    publishDir "$results_path/ngsplot"

    input:
    file input_files from bam_dedup_ngsplotconfig.toSortedList()

    output:
    file 'ngsplot_config' into ngsplot_config

    """
    #!/usr/bin/env Rscript
    datafiles = c( "${(input_files as List).join('", "')}" )
    table<-matrix(0,nrow=length(datafiles),ncol=3)
    table<-as.data.frame(table)
    for (i in 1:length(datafiles)){
       table[i,1]<-datafiles[i]
       table[i,2]<-(-1)
       tmp='\"'
       table[i,3]<-paste(tmp,gsub(".dedup.sort.bam.*\$","",as.character(datafiles[i])),tmp,sep="")
    }
    write.table(table,file="ngsplot_config",sep='\\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
    """
}



/*
* STEP 8 Ngsplot
*/

process ngsplot {

    errorStrategy 'ignore'

    module 'bioinfo-tools'
    module 'samtools'
    module 'R/3.2.3'
    module 'ngsplot'

    cpus 2
    memory '16 GB'
    time '120h'

    publishDir "$results_path/ngsplot"

    input:
    file input_bam_files from bam_dedup_ngsplot.toSortedList()
    file input_bai_files from bai_dedup_ngsplot.toSortedList()
    file ngsplot_config from ngsplot_config

    def REF
    if (params.genome=='GRCh37'){
        REF='hg19'
    }
    else if (params.genome=='GRCm38'){
        REF='mm10'
    }
    else if (params.genome=='NCBIM37'){
        REF='mm9'
    }
    else{
        error "No reference available for ngsplot!"
    }

    output:
    file '*.pdf' into results
    file '*_TSS' into results
    file '*_Gene' into results
    file '*.cnt' into results

    """
    /sw/apps/bioinfo/ngsplot/2.61/milou/bin/ngs.plot.r -G ${REF} -R genebody -C ${ngsplot_config} -O Genebody -D ensembl -FL 300

    /sw/apps/bioinfo/ngsplot/2.61/milou/bin/ngs.plot.r -G ${REF} -R tss -C ${ngsplot_config} -O TSS -FL 300
    """
}

/*
* STEP 9 MACS
*/

process macs {

    module 'bioinfo-tools'
    module 'MACS'
    module 'samtools'

    cpus 2
    memory '16 GB'
    time '24h'

    publishDir "$results_path/macs"

    input:
    file bam_for_macs from bam_dedup_macs.toSortedList()
    file bai_for_macs from bai_dedup_macs.toSortedList()

    set ChIPSampleID, CtrlSampleID, AnalysisID from MACSpara

    def REF
    if (params.genome=='GRCh37'){
        REF='hs'
    }
    else if (params.genome=='GRCm38'){
        REF='mm'
    }
    else if (params.genome=='NCBIM37'){
        REF='mm'
    }
    else{
        error "No reference available for MACS!"
    }


    output:
    file '*.bed' into results
    file '*.xls' into results
    file '*.r' into results
    file '*.narrowPeak' into results

    script:
    if (CtrlSampleID=='')
    """
        macs2 callpeak -t ${ChIPSampleID}.dedup.sort.bam -f BAM -g ${REF} -n ${AnalysisID} -q 0.01
    """

    else
    """
        macs2 callpeak -t ${ChIPSampleID}.dedup.sort.bam -c ${CtrlSampleID}.dedup.sort.bam -f BAM -g ${REF} -n ${AnalysisID} -q 0.01
    """
}

/*
 * STEP 10 MultiQC
 */

process multiqc {
    module 'bioinfo-tools'
    module 'MultiQC'

    memory '4GB'
    time '4h'

    publishDir "$results_path/MultiQC"

    errorStrategy 'ignore'

    input:
    file 'fastqc_report' from fastqc_html
    file 'trim_galore_report' from trim_galore_report
    file 'picard_report' from picard_report

    output:
    file 'multiqc_report.html'

     """
    multiqc -f  $PWD/results
    """
}






/*
 * Helper function, given a file Path
 * returns the file name region matching a specified glob pattern
 * starting from the beginning of the name up to last matching group.
 *
 * For example:
 *   readPrefix('/some/data/file_alpha_1.fa', 'file*_1.fa' )
 *
 * Returns:
 *   'file_alpha'
 */

def readPrefix( Path actual, template ) {

    final fileName = actual.getFileName().toString()

    def filePattern = template.toString()
    int p = filePattern.lastIndexOf('/')
    if( p != -1 ) filePattern = filePattern.substring(p+1)
    if( !filePattern.contains('*') && !filePattern.contains('?') )
        filePattern = '*' + filePattern

    def regex = filePattern
                    .replace('.','\\.')
                    .replace('*','(.*)')
                    .replace('?','(.?)')
                    .replace('{','(?:')
                    .replace('}',')')
                    .replace(',','|')

    def matcher = (fileName =~ /$regex/)
    if( matcher.matches() ) {
        def end = matcher.end(matcher.groupCount() )
        def prefix = fileName.substring(0,end)
        while(prefix.endsWith('-') || prefix.endsWith('_') || prefix.endsWith('.') )
          prefix=prefix[0..-2]

        return prefix
    }
    println(fileName)
    return fileName
}
