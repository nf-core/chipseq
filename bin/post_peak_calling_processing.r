#!/usr/bin/env Rscript

# Command line arguments
args = commandArgs(trailingOnly=TRUE)

R_lib <- as.character(args[1])
ref <- as.character(args[2])
Blacklist <- as.character(args[3])
input <- as.character(args[4:length(args)])

# Load / install required packages
.libPaths( c( .libPaths(), R_lib ) )

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

if (!require("biomaRt")){
    source("http://bioconductor.org/biocLite.R")
    biocLite("biomaRt", suppressUpdates=TRUE)
    library("biomaRt")
}

# Load / install required annotation databases
if(ref=="hs"){
	if (!require("BSgenome.Hsapiens.UCSC.hg19")){
        source("http://bioconductor.org/biocLite.R")
        biocLite("BSgenome.Hsapiens.UCSC.hg19", suppressUpdates=TRUE)
        library("BSgenome.Hsapiens.UCSC.hg19")
    }
    if (!require("org.Hs.eg.db")){
        source("http://bioconductor.org/biocLite.R")
        biocLite("org.Hs.eg.db", suppressUpdates=TRUE)
        library("org.Hs.eg.db")
    }
	data(TSS.human.GRCh37)
	annoData <- annoGR(TSS.human.GRCh37, feature="gene")
	orgAnnData <-"org.Hs.eg.db"
	mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
} else if(ref=="mm"){
	if (!require("BSgenome.Mmusculus.UCSC.mm10")){
        source("http://bioconductor.org/biocLite.R")
        biocLite("BSgenome.Mmusculus.UCSC.mm10", suppressUpdates=TRUE)
        library("BSgenome.Mmusculus.UCSC.mm10")
    }
    if (!require("org.Mm.eg.db")){
        source("http://bioconductor.org/biocLite.R")
        biocLite("org.Mm.eg.db", suppressUpdates=TRUE)
        library("org.Mm.eg.db")
    }
	data(TSS.mouse.GRCm38)
	annoData <- annoGR(TSS.mouse.GRCm38, feature="gene")
	orgAnnData <-"org.Mm.eg.db"
	mart=useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
} else{
	stop("Genome not supported!", call.=FALSE)
}

# Read in blacklist file and convert into range object
if (Blacklist!="No-filtering") {
	blacklist<-read.table(Blacklist,header=FALSE)
    blacklist_range<-with(blacklist,GRanges(sub("chr","",V1),IRanges(start=V2,end=V3),strand=Rle(rep("+",nrow(blacklist)))))
}

# Process output files from MACS: filtering peaks that overlap with blacklisted regions and annotating peaks
for (i in 1:length(input)) {
	
    # Read in raw peaks from MACS and convert into range object
    data<-read.table(input[i],header=TRUE)
    data_range<-with(data,GRanges(chr,IRanges(start=start,end=end),strand=Rle(rep("+",nrow(data))),length=length,summit=abs_summit,pileup=pileup,pvalue=X.log10.pvalue.,fold_enrichment=fold_enrichment,qvalue=X.log10.qvalue.,id=name))

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

    # Adding gene symbol
    if(class(try(final_anno<-addGeneIDs(annotatedPeak=final_anno,orgAnn=orgAnnData,IDs2Add="symbol")))=="try-error"){
	    feature_ids <- unique(final_anno$feature)
	    feature_ids <- feature_ids[!is.na(feature_ids)]
	    feature_ids <- feature_ids[feature_ids!=""]
	
	    IDs2Add<-getBM(attributes=c("ensembl_gene_id","external_gene_name"), filters = "ensembl_gene_id", values = feature_ids, mart=mart)
	
	    duplicated_ids<-IDs2Add[duplicated(IDs2Add[,"ensembl_gene_id"]),"ensembl_gene_id"]

	    if(length(duplicated_ids)>0){
	        IDs2Add.duplicated<-IDs2Add[IDs2Add[,"ensembl_gene_id"] %in% duplicated_ids,]
	        IDs2Add.duplicated<-condenseMatrixByColnames(as.matrix(IDs2Add.duplicated),"ensembl_gene_id")
	        IDs2Add<-IDs2Add[!(IDs2Add[,"ensembl_gene_id"] %in% duplicated_ids),]
	        IDs2Add<-rbind(IDs2Add,IDs2Add.duplicated)
	    }
	
	    final_anno$external_gene_name<-IDs2Add[match(IDs2Add$ensembl_gene_id,final_anno$feature),]$external_gene_name
    }
    else{
	    final_anno<-addGeneIDs(annotatedPeak=final_anno,orgAnn=orgAnnData,IDs2Add="symbol")
    }

    # Write annotated peaks to file
    final_anno_df<-as.data.frame(final_anno)
    newfilename<-paste(sub("_peaks.xls","",basename(input[i])),filter_flag,"_annotated.txt",sep="")
    write.table(final_anno_df,file=newfilename,quote=FALSE,sep="\t",eol="\n")
}
