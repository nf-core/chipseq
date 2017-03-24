#!/usr/bin/env Rscript

# Command line argument processing
datafiles <- commandArgs(trailingOnly=TRUE)

# Check input args
if (length(datafiles) < 1) {
  stop("Usage: ngs_config_generate.r [input-1] [input-2]..", call.=FALSE)
}

table<-matrix(0, nrow=length(datafiles), ncol=3)
table<-as.data.frame(table)
for (i in 1:length(datafiles)){
   table[i,1]<-datafiles[i]
   table[i,2]<-(-1)
   table[i,3]<-paste('"', gsub(".dedup.sorted.bam.*$", "", as.character(datafiles[i])), '"', sep="")
}
table<-table[order(table[,1]),]
write.table(table, file="ngsplot_config",sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
