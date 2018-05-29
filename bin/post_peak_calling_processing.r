#!/usr/bin/env Rscript

# R scripts for processing MACS output files (.xls)
# Author @chuan-wang https://github.com/chuan-wang

# Command line arguments
args <- commandArgs(trailingOnly=TRUE)

R_lib <- as.character(args[1])
ref <- as.character(args[2])
Blacklist <- as.character(args[3])
GTF <- as.character(args[4])
input <- as.character(args[5:length(args)])

# Load / install required packages
.libPaths( c( R_lib, .libPaths() ) )

if (!require("GenomicRanges")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges", suppressUpdates=TRUE)
  library("GenomicRanges")
}

if (!require("ChIPpeakAnno")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("ChIPpeakAnno", suppressUpdates=TRUE)
  library("ChIPpeakAnno")
}

if (!require("rtracklayer")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("rtracklayer", suppressUpdates=TRUE)
  library("rtracklayer")
}

if (!require("doParallel")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("doParallel", suppressUpdates=TRUE)
  library("doParallel")
}
if (!require("parallel")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("parallel", suppressUpdates=TRUE)
  library("parallel")
}

# Process annotation file
gtf <- import(GTF)

cores <- detectCores()/2
registerDoParallel(cores)

annotation <- as.data.frame(gtf)
annotation <- annotation[!duplicated(annotation),]
annotation <- makeGRangesFromDataFrame(annotation,keep.extra.columns = T)

anno <- plyr::adply(unique(annotation$gene_id),1, function(x){
  aux <- reduce(annotation[annotation$gene_id == x])
  aux$symbol <- unique(annotation[annotation$gene_id == x]$gene_name)
  as.data.frame(aux)
}, .id = NULL, .parallel = TRUE)

annoData <- makeGRangesFromDataFrame(anno, keep.extra.columns = T)

# Read in blacklist file and convert into range object
if (Blacklist!="No-filtering") {
  blacklist <- read.table(Blacklist,header=FALSE)
  blacklist_range <- with(blacklist,GRanges(sub("chr","",V1),IRanges(start=V2,end=V3),strand=Rle(rep("+",nrow(blacklist)))))
}

# Process output files from MACS: filtering peaks that overlap with blacklisted regions and annotating peaks
for (i in 1:length(input)) {
  
  # Avoid error that result file with 0 peak identified
  if(class(try(read.table(input[i],header=TRUE),silent=TRUE))=="try-error"){
    next
  }
  else{
    # Read in raw peaks from MACS and convert into range object
    data<-read.table(input[i],header=TRUE)
    data_range<-with(data,GRanges(chr,IRanges(start=start,end=end),strand=Rle(rep("+",nrow(data))),length=length,pileup=pileup,pvalue=X.log10.pvalue.,fold_enrichment=fold_enrichment,qvalue=X.log10.qvalue.,id=name))
    
    # Filtering peaks that overlap with blacklisted regions
    if (Blacklist!="No-filtering"){
      final<-data_range[data_range %outside% blacklist_range]
      filter_flag<-"_filtered"
    } else{
      final<-data_range
      filter_flag<-""
    }
    
    # Write peaks to txt and bed files
    final_df<-as.data.frame(final)
    newfilename<-paste(sub("_peaks.xls","",basename(input[i])),filter_flag,".txt",sep="")
    write.table(final_df,file=newfilename,quote=FALSE,sep="\t",eol="\n")
    
    df<-data.frame(seqnames=seqnames(final),starts=start(final)-1,ends=end(final),names=c(rep(".",length(final))),scores=c(rep(".",length(final))),strands=strand(final))
    newfilename<-paste(sub("_peaks.xls","",basename(input[i])),filter_flag,".bed",sep="")
    write.table(df, file=newfilename, quote=F, sep="\t", row.names=F, col.names=F)
    
    # Annotation
    final_anno<-annotatePeakInBatch(final, AnnotationData=annoData,output="overlapping", maxgap=5000L)
    
    # Write annotated peaks to file
    final_anno_df<-as.data.frame(final_anno)
    final_anno_df$gene_symbol<-anno[match(final_anno_df$feature,anno$gene),]$symbol
    newfilename<-paste(sub("_peaks.xls","",basename(input[i])),filter_flag,"_annotated.txt",sep="")
    write.table(final_anno_df,file=newfilename,quote=FALSE,sep="\t",eol="\n")
  }
}

# Show software versions
sessionInfo()