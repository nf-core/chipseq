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
data$peak<-""
data[data$x==crosscorr$phantom.cc$x,]$peak="phantom"
data[data$x==cc.peak[cc.peak$y==max(cc.peak$y),]$x,]$peak="peak"

write.csv(data, file=args[2], quote=FALSE, row.names=FALSE)