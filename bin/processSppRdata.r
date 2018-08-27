#!/usr/bin/env Rscript

# R scripts for convert the Rdata output file from phantompeakqualtools for strand-shift plot
# Author @chuan-wang https://github.com/chuan-wang

# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)

# Check input args
if (length(args) != 2) {
  stop("Usage: processSppRdata.r [ input.Rdata ] [ output.csv ]", call.=FALSE)
}

load(args[1])
data<-crosscorr$cross.correlation

write.table(data, file=args[2], sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE)