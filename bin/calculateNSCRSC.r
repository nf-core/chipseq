#!/usr/bin/env Rscript

# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)

# Check input args
if (length(args) != 1) {
  stop("Usage: calculateNSCRSC.r [ cross-correlation-file ]", call.=FALSE)
}

data<-read.table(args[1], header=FALSE)

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

write.table(data, file="cross_correlation_processed.txt", quote=FALSE, sep='\t', row.names=FALSE)