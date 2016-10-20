#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                    C H I P - S E Q   B E S T   P R A C T I C E
========================================================================================
 ChIP-seq Analysis Pipeline. Started May 2016.
 @Authors
 Chuan Wang <chuan.wang@scilifelab.se>
 Phil Ewels <phil.ewels@scilifelab.se>
 Release 2016-09-23
----------------------------------------------------------------------------------------
 Basic command:
 $ nextflow main.nf

 Pipeline variables can be configured with the following command line options:
 --genome [ID]
 --index [path to bwa index] (alternative to --genome)
 --reads [path to input files]
 --macsconfig [path to config file for MACS: line format: chip_sample_id,ctrl_sample_id,analysis_id]

 For example:
 $ nextflow main.nf -c ~/.nextflow/config --reads '*.fastq.gz' --macsconfig 'macssetup.config'
 $ nextflow main.nf -c ~/.nextflow/config --reads '*.R{1,2}.fastq.gz' --macsconfig 'macssetup.config'
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
version = 1.2

// Configurable variables
params.genome = 'GRCh37'
params.index = params.genomes[ params.genome ].bwa
params.name = "NGI ChIP-seq Best Practice"
params.reads = "data/*{1,2}*.fastq.gz"
params.macsconfig = "data/macsconfig"
params.extendReadsLen = 100
params.outdir = './results'

// R library locations
params.rlocation = "$HOME/R/nxtflow_libs/"
nxtflow_libs = file(params.rlocation)
nxtflow_libs.mkdirs()

single = 'null'

log.info "===================================="
log.info " ChIP-seq: v${version}"
log.info "===================================="
log.info "Reads        : ${params.reads}"
log.info "Genome       : ${params.genome}"
log.info "BWA Index    : ${params.index}"
log.info "MACS Config  : ${params.macsconfig}"
log.info "Extend Reads : ${params.extendReadsLen} bp"
log.info "Current home : $HOME"
log.info "Current user : $USER"
log.info "Current path : $PWD"
log.info "R libraries  : ${params.rlocation}"
log.info "Script dir   : $baseDir"
log.info "Working dir  : $workDir"
log.info "Output dir   : ${params.outdir}"
log.info "===================================="

// Validate inputs
index = file(params.index)
macsconfig = file(params.macsconfig)
if( !index.exists() ) exit 1, "Missing BWA index: '$index'"
if( !macsconfig.exists() ) exit 1, "Missing MACS config: '$macsconfig'. Specify path with --macsconfig"

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

read_files.into  { raw_reads_fastqc; raw_reads_trimgalore }


/*
 * Create a channel for macs config file
 */

macs_para = Channel
    .from(macsconfig.readLines())
    .map { line ->
        list = line.split(',')
        chip_sample_id = list[0]
        ctrl_sample_id = list[1]
        analysis_id = list[2]
        [ chip_sample_id, ctrl_sample_id, analysis_id ]
    }



/*
 * STEP 1 - FastQC
 */

process fastqc {

    tag "$name"

    module 'bioinfo-tools'
    module 'FastQC'

    memory { 2.GB * task.attempt }
    time { 4.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'warning' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    set val(name), file(reads:'*') from raw_reads_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc -q $reads
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

    cpus 2
    memory { 4.GB * task.attempt }
    time { 8.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/trimgalore", mode: 'copy'

    input:
    set val(name), file(reads:'*') from raw_reads_trimgalore

    output:
    file '*.gz' into trimmed_reads
    file '*trimming_report.txt' into trimgalore_results

    script:
    single = reads instanceof Path
    if(single) {
        """
        trim_galore --gzip $reads
        """
    } else {
        """
        trim_galore --paired --gzip $reads
        """
    }
}


/*
 * STEP 3.1 - align with bwa
 */

process bwa {

    tag "$reads"

    module 'bioinfo-tools'
    module 'bwa'
    module 'samtools'

    cpus 8
    memory { 64.GB * task.attempt }
    time { 48.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bwa", mode: 'copy'

    input:
    file reads from trimmed_reads

    output:
    file '*.bam' into bwa_bam
    stdout into bwa_logs

    script:
    """
    f='$reads';f=(\$f);f=\${f[0]};f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%.*1_val_1};f=\${f%_1_val_1};f=\${f%_trimmed}

    bwa mem -M $index $reads | samtools view -bT $index - > \${f}.bam
    """
}


/*
 * STEP 3.2 - post-alignment processing
 */

process samtools {

    tag "$bwa_bam"

    module 'bioinfo-tools'
    module 'samtools'
    module 'BEDTools'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 120.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bwa", mode: 'copy'

    input:
    file bwa_bam

    output:
    file '*.sorted.bam' into bam_picard
    file '*.sorted.bam.bai' into bwa_bai
    file '*.sorted.bed' into bed_total

    script:
    """
    f='$bwa_bam';f=\${f%.bam}

    samtools sort \${f}.bam \${f}.sorted
    samtools index \${f}.sorted.bam
    bedtools bamtobed -i \${f}.sorted.bam | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 > \${f}.sorted.bed
    """
}


/*
 * STEP 4 Picard
 */

process picard {

    tag "$bam_picard"

    module 'bioinfo-tools'
    module 'picard/2.0.1'
    module 'BEDTools'
    module 'samtools'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 24.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/picard", mode: 'copy'

    input:
    file bam_picard

    output:
    file '*.dedup.sorted.bam' into bam_dedup_spp, bam_dedup_ngsplotconfig, bam_dedup_ngsplot, bam_dedup_deepTools, bam_dedup_macs
    file '*.dedup.sorted.bam.bai' into bai_dedup_deepTools, bai_dedup_ngsplot, bai_dedup_macs
    file '*.dedup.sorted.bed' into bed_dedup
    file '*.picardDupMetrics.txt' into picard_reports

    script:
    """
    f='$bam_picard';f=\${f%.sorted.bam}

    java -Xmx2g -jar \$PICARD_HOME/picard.jar MarkDuplicates \\
        INPUT=$bam_picard \\
        OUTPUT=\${f}.dedup.bam \\
        ASSUME_SORTED=true \\
        REMOVE_DUPLICATES=true \\
        METRICS_FILE=\${f}.picardDupMetrics.txt \\
        VALIDATION_STRINGENCY=LENIENT \\
        PROGRAM_RECORD_ID='null'

    samtools sort \${f}.dedup.bam \${f}.dedup.sorted
    samtools index \${f}.dedup.sorted.bam
    bedtools bamtobed -i \${f}.dedup.sorted.bam | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 > \${f}.dedup.sorted.bed
    """
}


/*
 * STEP 5 Read_count_statistics
 */

process countstat {

    cpus 2
    memory { 4.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/countstat", mode: 'copy'

    input:
    file input from bed_total .mix(bed_dedup) .toSortedList()

    output:
    file 'read_count_statistics.txt' into countstat_results

    """
    #! /usr/bin/env perl

    use strict;
    use warnings;

    open(OUTPUT, ">read_count_statistics.txt");
    my @fileList = qw($input);

    print OUTPUT "File\\tTotalCounts\\tUniqueCounts\\tUniqueStartCounts\\tUniqueRatio\\tUniqueStartRatio\\n";

    foreach my \$f(@fileList){

     open(IN,"<\$f")||die \$!;
     my \$Tcnt=0;
     my \$prev="NA";
     my \$lcnt=0;
     my \$Tcnt_2=0;
     my \$prev_2="NA";

     while(<IN>){
       chomp;
       my @line=split("\\t",\$_);
       \$lcnt++;
       my \$t = join("_",@line[0..2]);
       \$Tcnt++ unless(\$t eq \$prev);
       \$prev=\$t;

       my \$t_2 = join("_",@line[0..1]);
       \$Tcnt_2++ unless(\$t_2 eq \$prev_2);
       \$prev_2=\$t_2;
     }

     print OUTPUT "\$f\\t\$lcnt\\t\$Tcnt\\t\$Tcnt_2\\t".(\$Tcnt/\$lcnt)."\\t".(\$Tcnt_2/\$lcnt)."\\n";
     close(IN);
    }
    close(OUTPUT);
    """
}


/*
 * STEP 6.1 Phantompeakqualtools
 */

process phantompeakqualtools {

    tag "$bam_dedup_spp"

    module 'bioinfo-tools'
    module 'samtools'
    module 'R/3.1.0'
    module 'phantompeakqualtools'
    module 'SPP'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 24.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/phantompeakqualtools", mode: 'copy'

    input:
    file bam_dedup_spp

    output:
    file '*.pdf' into spp_results
    file '*.spp.out' into spp_out, spp_out_mqc

    script:
    """
    f='$bam_dedup_spp';f=\${f%.dedup.sorted.bam}
    run_spp.R -c=$bam_dedup_spp -savp -out=\${f}.spp.out
    """
}


/*
 * STEP 6.2 Combine_spp_out
 */

process combinesppout {

    cpus 1
    memory { 2.GB * task.attempt }
    time { 1.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/phantompeakqualtools", mode: 'copy'

    input:
    file spp_out_list from spp_out.toSortedList()

    output:
    file 'cross_correlation.txt' into cross_correlation

    script:
    """
    cat $spp_out_list > cross_correlation.txt
    """
}


/*
 * STEP 6.3 Calculate_NSC_RSC
 */

process calculateNSCRSC {

    cpus 1
    memory { 2.GB * task.attempt }
    time { 1.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/phantompeakqualtools", mode: 'copy'

    input:
    file cross_correlation

    output:
    file 'crosscorrelation_processed.txt' into calculateNSCRSC_results

    script:
    """
    #!/usr/bin/env Rscript

    data<-read.table("${cross_correlation}",header=FALSE)

    data[,12]<-NA
    data[,13]<-NA
    data[,14]<-NA
    data[,15]<-NA

    colnames(data)[14]<-"NSC"
    colnames(data)[15]<-"RSC"

    for (i in 1:nrow(data)){
	       data[i,12]<-as.numeric(unlist(strsplit(as.character(data[i,4]),","))[1])
	       data[i,13]<-as.numeric(unlist(strsplit(as.character(data[i,6]),","))[1])
	       data[i,14]<-round(data[i,12]/as.numeric(data[i,8]),2)
	       data[i,15]<-round((data[i,12]-as.numeric(data[i,8]))/(data[i,13]-as.numeric(data[i,8])),2)
    }

    write.table(data, file="crosscorrelation_processed.txt", quote=FALSE, sep='\\t', row.names=FALSE)
    """
}


/*
 * STEP 7 deepTools
 */

process deepTools {

    module 'bioinfo-tools'
    module 'deepTools'

    cpus 4
    memory { 32.GB * task.attempt }
    time { 120.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/deepTools", mode: 'copy'

    input:
    file bam from bam_dedup_deepTools.toSortedList()
    file bai from bai_dedup_deepTools.toSortedList()

    output:
    file 'multiBamSummary.npz' into deepTools_bamsummary
    file '*.{pdf,png}' into deepTools_results

    script:
    """
    plotFingerprint \\
        -b $bam \\
        --plotFile fingerprints.pdf \\
        --extendReads=${params.extendReadsLen} \\
        --skipZeros \\
        --ignoreDuplicates \\
        --numberOfSamples 50000 \\
        --binSize=500 \\
        --plotFileFormat=pdf \\
        --plotTitle="Fingerprints"

    multiBamSummary \\
        bins \\
        --binSize=10000 \\
        --bamfiles $bam \\
        -out multiBamSummary.npz \\
        --extendReads=${params.extendReadsLen} \\
        --ignoreDuplicates \\
        --centerReads

    plotCorrelation \\
        -in multiBamSummary.npz \\
        -o scatterplot_PearsonCorr_multiBamSummary.png \\
        --corMethod pearson \\
        --skipZeros \\
        --removeOutliers \\
        --plotTitle "Pearson Correlation of Read Counts" \\
        --whatToPlot scatterplot \\

    plotCorrelation \\
        -in multiBamSummary.npz \\
        -o heatmap_SpearmanCorr_multiBamSummary.png \\
        --corMethod spearman \\
        --skipZeros \\
        --plotTitle "Spearman Correlation of Read Counts" \\
        --whatToPlot heatmap \\
        --colorMap RdYlBu \\
        --plotNumbers \\
    """
}


/*
 * STEP 8.1 Generate config file for ngsplot
 */

process ngs_config_generate {

    memory { 2.GB * task.attempt }
    time { 1.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/ngsplot", mode: 'copy'

    input:
    file input_files from bam_dedup_ngsplotconfig.toSortedList()

    output:
    file 'ngsplot_config' into ngsplot_config

    script:
    """
    #!/usr/bin/env Rscript

    datafiles = c( "${(input_files as List).join('", "')}" )
    table<-matrix(0, nrow=length(datafiles), ncol=3)
    table<-as.data.frame(table)
    for (i in 1:length(datafiles)){
       table[i,1]<-datafiles[i]
       table[i,2]<-(-1)
       tmp='\"'
       table[i,3]<-paste(tmp, gsub(".dedup.sorted.bam.*\$", "", as.character(datafiles[i])), tmp, sep="")
    }
    table<-table[order(table[,1]),]
    write.table(table, file="ngsplot_config",sep='\\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
    """
}


/*
 * STEP 8.2 Ngsplot
 */

process ngsplot {

    module 'bioinfo-tools'
    module 'samtools'
    module 'R/3.2.3'
    module 'ngsplot'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 120.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/ngsplot", mode: 'copy'

    input:
    file input_bam_files from bam_dedup_ngsplot.toSortedList()
    file input_bai_files from bai_dedup_ngsplot.toSortedList()
    file ngsplot_config from ngsplot_config

    output:
    file '*.{pdf}' into ngsplot_results

    script:

    def REF
    if (params.genome == 'GRCh37'){ REF = 'hg19' }
    else if (params.genome == 'GRCm38'){ REF = 'mm10' }
    else { error "No reference / reference not supported available for ngsplot! >${params.genome}<" }

    """
    \$NGSPLOT/bin/ngs.plot.r \\
        -G $REF \\
        -R genebody \\
        -C $ngsplot_config \\
        -O Genebody \\
        -D ensembl \\
        -FL 300

    \$NGSPLOT/bin/ngs.plot.r \\
        -G $REF \\
        -R tss \\
        -C $ngsplot_config \\
        -O TSS \\
        -FL 300
    """
}


/*
 * STEP 9 MACS
 */

process macs {

    module 'bioinfo-tools'
    module 'MACS'
    module 'samtools'

    cpus 4
    memory { 32.GB * task.attempt }
    time { 24.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/macs", mode: 'copy'

    input:
    file bam_for_macs from bam_dedup_macs.toSortedList()
    file bai_for_macs from bai_dedup_macs.toSortedList()
    set chip_sample_id, ctrl_sample_id, analysis_id from macs_para

    output:
    file '*.{bed,xls,r,narrowPeak}' into macs_results

    script:

    def REF
    if (params.genome == 'GRCh37'){ REF = 'hs' }
    else if (params.genome == 'GRCm38'){ REF = 'mm' }
    else { error "No reference / reference not supported available for MACS! >${params.genome}<" }

    if (ctrl_sample_id == '') {
        ctrl = ''
    } else {
        ctrl = "-c ${ctrl_sample_id}.dedup.sorted.bam"
    }
    """
    macs2 callpeak \\
        -t ${chip_sample_id}.dedup.sorted.bam \\
        $ctrl \\
        -f BAM \\
        -g $REF \\
        -n $analysis_id \\
        -q 0.01
    """
}


/*
 * STEP 10 MultiQC
 */

process multiqc {
    module 'bioinfo-tools'
    module 'MultiQC'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 120.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file ('fastqc/*') from fastqc_results.toList()
    file ('trimgalore/*') from trimgalore_results.toList()
    file ('bwa/*') from bwa_logs.toList()
    file ('picard/*') from picard_reports.toList()
    file ('phantompeakqualtools/*') from spp_out_mqc.toList()

    output:
    file '*multiqc_report.html' into multiqc_report
    file '*multiqc_data' into multiqc_data

    """
    multiqc -f .
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
