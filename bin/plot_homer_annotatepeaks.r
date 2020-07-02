#!/usr/bin/env Rscript

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(ggplot2)
library(reshape2)
library(scales)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(make_option(c("-i", "--homer_files"), type="character", default=NULL, help="Comma-separated list of homer annotated text files.", metavar="path"),
                    make_option(c("-s", "--sample_ids"), type="character", default=NULL, help="Comma-separated list of sample ids associated with homer annotated text files. Must be unique and in same order as homer files input.", metavar="string"),
                    make_option(c("-o", "--outdir"), type="character", default='./', help="Output directory", metavar="path"),
                    make_option(c("-p", "--outprefix"), type="character", default='homer_annotation', help="Output prefix", metavar="string"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$homer_files)){
    print_help(opt_parser)
    stop("At least one homer annotated file must be supplied", call.=FALSE)
}
if (is.null(opt$sample_ids)){
    print_help(opt_parser)
    stop("Please provide sample ids associated with homer files.", call.=FALSE)
}

if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir,recursive=TRUE)
}

HomerFiles <- unlist(strsplit(opt$homer_files,","))
SampleIDs <- unlist(strsplit(opt$sample_ids,","))
if (length(HomerFiles) != length(SampleIDs)) {
    print_help(opt_parser)
    stop("Number of sample ids must equal number of homer annotated files.", call.=FALSE)
}

################################################
################################################
## READ IN DATA                               ##
################################################
################################################

plot.dat <- data.frame()
plot.dist.dat <- data.frame()
plot.feature.dat <- data.frame()
for (idx in 1:length(HomerFiles)) {

    sampleid = SampleIDs[idx]
    anno.dat <- read.csv(HomerFiles[idx], sep="\t", header=TRUE)
    anno.dat <- anno.dat[,c("Annotation","Distance.to.TSS","Nearest.PromoterID")]

    ## REPLACE UNASSIGNED FEATURE ENTRIES WITH SENSIBLE VALUES
    unassigned <- which(is.na(as.character(anno.dat$Distance.to.TSS)))
    anno.dat$Distance.to.TSS[unassigned] <- 1000000

    anno.dat$Annotation <- as.character(anno.dat$Annotation)
    anno.dat$Annotation[unassigned] <- "Unassigned"
    anno.dat$Annotation <- as.factor(anno.dat$Annotation)

    anno.dat$Nearest.PromoterID <- as.character(anno.dat$Nearest.PromoterID)
    anno.dat$Nearest.PromoterID[unassigned] <- "Unassigned"
    anno.dat$Nearest.PromoterID <- as.factor(anno.dat$Nearest.PromoterID)

    anno.dat$name <- rep(sampleid,nrow(anno.dat))
    anno.dat$Distance.to.TSS <- abs(anno.dat$Distance.to.TSS) + 1
    plot.dat <- rbind(plot.dat,anno.dat)

    ## GET ANNOTATION COUNTS
    anno.freq <- as.character(lapply(strsplit(as.character(anno.dat$Annotation)," "), function(x) x[1]))
    anno.freq <- as.data.frame(table(anno.freq))
    colnames(anno.freq) <- c("feature",sampleid)
    anno.melt <- melt(anno.freq)
    plot.feature.dat <- rbind(plot.feature.dat,anno.melt)

    ## GET CLOSEST INSTANCE OF GENE TO ANY GIVEN PEAK
    unique.gene.dat <- anno.dat[order(anno.dat$Distance.to.TSS),]
    unique.gene.dat <- unique.gene.dat[!duplicated(unique.gene.dat$Nearest.PromoterID), ]
    dist.freq <- rep("> 10kb",nrow(unique.gene.dat))
    dist.freq[which(unique.gene.dat$Distance.to.TSS < 10000)] <- "< 10kb"
    dist.freq[which(unique.gene.dat$Distance.to.TSS < 5000)] <- "< 5kb"
    dist.freq[which(unique.gene.dat$Distance.to.TSS < 2000)] <- "< 2kb"
    dist.freq <- as.data.frame(table(dist.freq))
    colnames(dist.freq) <- c("distance",sampleid)
    dist.melt <- melt(dist.freq)
    plot.dist.dat <- rbind(plot.dist.dat,dist.melt)

}
plot.dat$name <- factor(plot.dat$name, levels=sort(unique(as.character(plot.dat$name))))
plot.dist.dat$variable <- factor(plot.dist.dat$variable, levels=sort(unique(as.character(plot.dist.dat$variable))))
plot.feature.dat$variable <- factor(plot.feature.dat$variable, levels=sort(unique(as.character(plot.feature.dat$variable))))

summary.dat <- dcast(plot.feature.dat, variable ~ feature, value.var="value")
colnames(summary.dat)[1] <- "sample"
write.table(summary.dat,file=file.path(opt$outdir,paste(opt$outprefix,".summary.txt",sep="")),sep="\t",row.names=F,col.names=T,quote=F)

################################################
################################################
## PLOTS                                      ##
################################################
################################################

PlotFile <- file.path(opt$outdir,paste(opt$outprefix,".plots.pdf",sep=""))
pdf(PlotFile,height=6,width=3*length(HomerFiles))

## FEATURE COUNT STACKED BARPLOT
plot  <- ggplot(plot.feature.dat, aes(x=variable, y=value, group=feature)) +
         geom_bar(stat="identity", position = "fill", aes(colour=feature,fill=feature), alpha = 0.3) +
         xlab("") +
         ylab("% Feature") +
         ggtitle("Peak Location Relative to Annotation") +
         scale_y_continuous(labels = percent_format()) +
         theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text.y = element_text(colour="black"),
               axis.text.x= element_text(colour="black",face="bold"),
               axis.line.x = element_line(size = 1, colour = "black", linetype = "solid"),
               axis.line.y = element_line(size = 1, colour = "black", linetype = "solid"))
print(plot)

## DISTANCE TO CLOSEST GENE ACROSS ALL PEAKS STACKED BARPLOT
plot  <- ggplot(plot.dist.dat, aes(x=variable, y=value, group=distance)) +
         geom_bar(stat="identity", position = "fill", aes(colour=distance,fill=distance), alpha = 0.3) +
         xlab("") +
         ylab("% Unique genes to closest peak") +
         ggtitle("Distance of Closest Peak to Gene") +
         scale_y_continuous(labels = percent_format()) +
         theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text.y = element_text(colour="black"),
               axis.text.x= element_text(colour="black",face="bold"),
               axis.line.x = element_line(size = 1, colour = "black", linetype = "solid"),
               axis.line.y = element_line(size = 1, colour = "black", linetype = "solid"))
print(plot)

## VIOLIN PLOT OF PEAK DISTANCE TO TSS
plot  <- ggplot(plot.dat, aes(x=name, y=Distance.to.TSS)) +
         geom_violin(aes(colour=name,fill=name), alpha = 0.3) +
         geom_boxplot(width=0.1) +
         xlab("") +
         ylab(expression(log[10]*" distance to TSS")) +
         ggtitle("Peak Distribution Relative to TSS") +
         scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
         theme(legend.position="none",
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text.y = element_text(colour="black"),
               axis.text.x= element_text(colour="black",face="bold"),
               axis.line.x = element_line(size = 1, colour = "black", linetype = "solid"),
               axis.line.y = element_line(size = 1, colour = "black", linetype = "solid"))
print(plot)
dev.off()

################################################
################################################
################################################
################################################
