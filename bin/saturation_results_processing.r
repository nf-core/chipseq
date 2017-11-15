#!/usr/bin/env Rscript

# R scripts for summarizing saturation analysis results using MACS output files (.xls) as input.
# Version 1.0
# Author @chuan-wang https://github.com/chuan-wang

# Command line arguments
args <- commandArgs(trailingOnly=TRUE)

R_lib <- as.character(args[1])
config <- as.character(args[2])
countstat <- as.character(args[3])
input <- as.character(args[4:length(args)])

# Load / install required packages
.libPaths( c( R_lib, .libPaths() ) )

if (!require("ggplot2")){
    install.packages("ggplot2", dependencies=TRUE, repos='http://cloud.r-project.org/')
    library("ggplot2")
}

# Read in MACS config file
macsconfig <- read.csv(config,header=F)

# Read in statistics of read counts
readcount <- read.table(countstat,header=T)
readcount$Sample <- gsub("\\..*$","",as.character(readcount$File))
readcount$Category <- NA
for (i in 1:nrow(readcount)) {
	if (gsub("^[^.]*","",as.character(readcount$File[i])) == ".dedup.sorted.bed") {
		readcount$Category[i] <- "Duplicates_removed"
	} else {
		readcount$Category[i] <- "Raw"
	}
}
readcount_duprm<-readcount[readcount$Category=="Duplicates_removed",]

# Read in output files from MACS
meta<-matrix(,nrow=length(input),ncol=2)
meta<-as.data.frame(meta)
colnames(meta)<-c("InputFile","Peaks")

for (i in 1:length(input)) {
	
	if(class(try(read.table(input[i],header=TRUE),silent=TRUE))=="try-error"){
		# Result file with 0 peak identified
		meta$Peaks[i] <- 0
		meta$InputFile[i] <- sub("_peaks.xls","",basename(input[i]))
	}
	else{
		# Result file with >0 peaks identified
		tmp<-read.table(input[i],header=TRUE)
		meta$Peaks[i] <- nrow(tmp)
		meta$InputFile[i] <- sub("_peaks.xls","",basename(input[i]))
	}
}

# Fill in AnalysisID, ChIPSampleID, and reads used for peak calling
meta$AnalysisID <- gsub("[.][0-9][.][0-9]*$","",as.character(meta$InputFile))
meta$ReadProportion <- as.numeric(sub(".*([0-9][.][0-9])$","\\1",as.character(meta$InputFile)))
meta$ChIPSampleID <- as.character(macsconfig[match(meta$AnalysisID,as.character(macsconfig$V3)),]$V1)
meta$ReadNumber <- meta$ReadProportion*as.numeric(as.character(readcount_duprm[match(meta$ChIPSampleID,as.character(readcount_duprm$Sample)),]$TotalCounts))

# Write summary to txt file
write.table(meta,file="Saturation_analysis_summary.txt",quote=FALSE,sep="\t",eol="\n")

# Plot with input read number as x axis and peak number as y axis
pdf("Saturation_analysis_plot.pdf",width=24,height=12)
ggplot(data=meta,aes(ReadNumber,Peaks,col=AnalysisID))+geom_line(lwd=2)+xlab("Input Reads")+ylab("Peaks Identified")+ggtitle("Saturation Plot")+ theme(plot.title = element_text(size=rel(2), face="bold"),axis.text.x=element_text(size=rel(1.5)),axis.text.y=element_text(size=rel(1.5)),axis.title.x=element_text(size=rel(1.5)),axis.title.y=element_text(size=rel(1.5)),legend.text=element_text(size=rel(1.5)),legend.title=element_blank())
dev.off()