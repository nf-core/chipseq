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

# Process annotation file
gtf<-import(GTF)
annotation<-as.data.frame(gtf)
annotation<-annotation[!duplicated(annotation),]
genes<-unique(as.character(annotation$gene_id))
anno<-matrix(nrow=length(genes),ncol=6)
anno<-as.data.frame(anno)
colnames(anno)<-c("chr","start","end","strand","gene","symbol")
anno$gene<-genes

for(i in 1:nrow(anno)){
  anno$chr[i]<-unique(as.character(annotation[as.character(annotation$gene_id)==anno$gene[i],]$seqnames))
  anno$start[i]<-min(as.numeric(as.character(annotation[as.character(annotation$gene_id)==anno$gene[i],]$start)))
  anno$end[i]<-max(as.numeric(as.character(annotation[as.character(annotation$gene_id)==anno$gene[i],]$end)))
  anno$strand[i]<-unique(as.character(annotation[as.character(annotation$gene_id)==anno$gene[i],]$strand))
  anno$symbol[i]<-unique(as.character(annotation[as.character(annotation$gene_id)==anno$gene[i],]$gene_name))
}

annoData<-with(anno,GRanges(seqnames=chr,ranges=IRanges(start=start,end=end,names=gene),strand=strand,symbol=symbol))

# Read in blacklist file and convert into range object
if (Blacklist!="No-filtering") {
	blacklist<-read.table(Blacklist,header=FALSE)
    blacklist_range<-with(blacklist,GRanges(sub("chr","",V1),IRanges(start=V2,end=V3),strand=Rle(rep("+",nrow(blacklist)))))
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