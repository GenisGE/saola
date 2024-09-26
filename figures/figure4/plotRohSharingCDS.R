source("/home/lpq293/mypopgen/saola/info/loadPopInfo.R")

classifyRohs <- function(nrohs, nhaps, length){
    
    # vectors with for each segment number of rohs, number of different haps, length of rohs
    length <- as.numeric(length)
    sapply(0:max(nrohs), function(x) c(sum(length[nrohs==x & nhaps==min(1,x)]),
                                       sum(length[nrohs==x & nhaps > 1])) / sum(length) )
    
}



f <- "/home/lpq293/mypopgen/saola/analyses/rohsharing/results/rohsharing_tables/rohsharing_summary_withcdsprop.tsv"

dat <- read.table(f, h=T, stringsAsFactors=F)

k <- !is.na(dat$nhaps)
lost <- sum(dat$length[!k])
kept <- sum(as.numeric(dat$length[k]))
dat <- dat[k,]


northerncols <- rev(c("#bcd1e2", "#385f82"))
southerncols <- rev(c("#ffe199",  "#b37d00"))
combinedcols <- rev(c("#cccccc", "#5d5d5d"))
saola_colors2 <- c("Northern" = "#325675", "Southern" = "#7f5900", "Combined" = "#8f8f8f")
saola_colors2 <- c(saola_colors, "Combined"= "#8f8f8f")


outpng <- "nicePlotRohSharingCDS2.png"
bitmap(outpng, w=4,h=1.5,res=300)
par(mfrow=c(2,3), mar=c(2.1, 0, 2.1, 0), oma=c(1,6,0,0))

all_a <- c(classifyRohs(dat$num, dat$nhaps, dat$length),
           classifyRohs(dat$nrohNorth, dat$nhapsNorth, dat$length),
           classifyRohs(dat$nrohSouth, dat$nhapsSouth, dat$length))

ylims <- c(0, max(all_a))
a <- classifyRohs(dat$num, dat$nhaps, dat$length)
bppos <- barplot(a, xaxt="n", ylab="Genome fraction",
        main="All samples", cex.lab=cexlab, cex.main=cexmain, ylim=ylims, xpd=NA,
        col=combinedcols, cex.axis=cexaxis)
title( xlab="Number of samples in ROH", cex.lab=cexlab2, line=1.75, xpd=NA)
axis(1, at=bppos, labels=0:(ncol(a)-1), tick=F, line=-0.5, cex.axis=cexaxis)
legend("topright", fill=combinedcols, 
                legend=c("Same haplotype", "Different haplotypes"), bty="n", 
                cex=1.25)

a <- classifyRohs(dat$nrohNorth, dat$nhapsNorth, dat$length)
bppos <- barplot(a, xaxt="n", ylab="",
        main="Northern samples", cex.lab=1.5, cex.main=cexmain, ylim=ylims, yaxt="n", xpd=NA,
        col=northerncols)
title( xlab="Number of samples in ROH", cex.lab=cexlab2, line=1.75, xpd=NA)
axis(1, at=bppos, labels=0:(ncol(a)-1), tick=F, line=-0.5, cex.axis=cexaxis)


a <- classifyRohs(dat$nrohSouth, dat$nhapsSouth, dat$length)
bppos <- barplot(a, xaxt="n", ylab="",
        main="Southern samples", cex.lab=1.5, cex.main=cexmain, ylim=ylims, yaxt="n", xpd=NA,
        col=southerncols)
title(xlab="Number of samples in ROH", cex.lab=cexlab2, line=1.75, xpd=NA)
axis(1, at=bppos, labels=0:(ncol(a)-1), tick=F, line=-0.5, cex.axis=cexaxis)


all_a <- max(unlist(sapply(list(dat$num, dat$nrohNorth, dat$nrohSout),
                function(y)  tapply(1:nrow(dat), y, function(x) sum(dat$length[x] * dat$cds_proportion[x]) / sum(dat$length[x])))))

ylims <- c(max(all_a), 0)

a <- tapply(1:nrow(dat), dat$num, function(x) sum(dat$length[x] * dat$cds_proportion[x]) / sum(dat$length[x]))

bppos <- barplot(a, names.arg=names(a), ylab="Proporiton of\ncoding sequence",
                 main= "", cex.lab=cexlab, cex.main=1.5, cex.names=1.2, xlab="", xaxt="n",
                 ylim=ylims, xpd=NA,
                 col=saola_colors2["Combined"], cex.axis=cexaxis)
axis(3, at=bppos, labels=names(a), tick=F, line=-1, cex.axis=cexaxis)
abline(h=sum(as.numeric(dat$length * dat$cds_proportion)) / sum(as.numeric(dat$length)), lwd=2, lty=2, col="darkred")
legend(x=3,y=0.0136, lwd=2, lty=2, 
        col="darkred", legend="Genome-wide proportion\nof coding sequence", 
        cex=1.4,
        bty="n", xpd=NA)

a <- tapply(1:nrow(dat), dat$nrohNorth, function(x) sum(dat$length[x] * dat$cds_proportion[x]) / sum(dat$length[x]))

bppos <- barplot(a, names.arg=names(a), ylab="",
                 main= "", cex.lab=1.5, cex.main=1.5, cex.names=1.2, xlab="", xaxt="n", yaxt="n",
                 ylim=ylims,
                 col=saola_colors2["Northern"])
axis(3, at=bppos, labels=names(a), tick=F, line=-1, cex.axis=cexaxis)
abline(h=sum(as.numeric(dat$length * dat$cds_proportion)) / sum(as.numeric(dat$length)), lwd=2, lty=2, col="darkred")

a <- tapply(1:nrow(dat), dat$nrohSouth, function(x) sum(dat$length[x] * dat$cds_proportion[x]) / sum(dat$length[x]))

bppos <- barplot(a, names.arg=names(a), ylab="",
                 main= "", cex.lab=1.5, cex.main=1.5, cex.names=1.2, xlab="", xaxt="n", yaxt="n",
                 ylim=ylims,
                 col=saola_colors2["Southern"])
axis(3, at=bppos, labels=names(a), tick=F, line=-1, cex.axis=cexaxis)
abline(h=sum(as.numeric(dat$length * dat$cds_proportion)) / sum(as.numeric(dat$length)), lwd=2, lty=2, col="darkred")
dev.off()




cexmain <- 1.7
cexlab <- 1.7
cexlab2 <- 1.65
cexaxis <- 1.4

outpdf <- "nicePlotRohSharingCDS2.pdf"
pdf(outpdf, w=4 * 2,h=1.5 * 2)
par(mfrow=c(2,3), mar=c(2.1, 0, 2.1, 0), oma=c(1,6,0,0))

all_a <- c(classifyRohs(dat$num, dat$nhaps, dat$length),
           classifyRohs(dat$nrohNorth, dat$nhapsNorth, dat$length),
           classifyRohs(dat$nrohSouth, dat$nhapsSouth, dat$length))

ylims <- c(0, max(all_a))
a <- classifyRohs(dat$num, dat$nhaps, dat$length)
bppos <- barplot(a, xaxt="n", ylab="Genome fraction",
        main="All samples", cex.lab=cexlab, cex.main=cexmain, ylim=ylims, xpd=NA,
        col=combinedcols, cex.axis=cexaxis)
title( xlab="Number of samples in ROH", cex.lab=cexlab2, line=1.75, xpd=NA)
axis(1, at=bppos, labels=0:(ncol(a)-1), tick=F, line=-0.5, cex.axis=cexaxis)
legend("topright", fill=combinedcols, 
                legend=c("Same haplotype", "Different haplotypes"), bty="n", 
                cex=1.25)

a <- classifyRohs(dat$nrohNorth, dat$nhapsNorth, dat$length)
bppos <- barplot(a, xaxt="n", ylab="",
        main="Northern samples", cex.lab=1.5, cex.main=cexmain, ylim=ylims, yaxt="n", xpd=NA,
        col=northerncols)
title( xlab="Number of samples in ROH", cex.lab=cexlab2, line=1.75, xpd=NA)
axis(1, at=bppos, labels=0:(ncol(a)-1), tick=F, line=-0.5, cex.axis=cexaxis)


a <- classifyRohs(dat$nrohSouth, dat$nhapsSouth, dat$length)
bppos <- barplot(a, xaxt="n", ylab="",
        main="Southern samples", cex.lab=1.5, cex.main=cexmain, ylim=ylims, yaxt="n", xpd=NA,
        col=southerncols)
title(xlab="Number of samples in ROH", cex.lab=cexlab2, line=1.75, xpd=NA)
axis(1, at=bppos, labels=0:(ncol(a)-1), tick=F, line=-0.5, cex.axis=cexaxis)


all_a <- max(unlist(sapply(list(dat$num, dat$nrohNorth, dat$nrohSout),
                function(y)  tapply(1:nrow(dat), y, function(x) sum(dat$length[x] * dat$cds_proportion[x]) / sum(dat$length[x])))))

ylims <- c(max(all_a), 0)

a <- tapply(1:nrow(dat), dat$num, function(x) sum(dat$length[x] * dat$cds_proportion[x]) / sum(dat$length[x]))

bppos <- barplot(a, names.arg=names(a), ylab="Proporiton of\ncoding sequence",
                 main= "", cex.lab=cexlab, cex.main=1.5, cex.names=1.2, xlab="", xaxt="n",
                 ylim=ylims, xpd=NA,
                 col=saola_colors2["Combined"], cex.axis=cexaxis)
axis(3, at=bppos, labels=names(a), tick=F, line=-1, cex.axis=cexaxis)
abline(h=sum(as.numeric(dat$length * dat$cds_proportion)) / sum(as.numeric(dat$length)), lwd=2, lty=2, col="darkred")
legend(x=3,y=0.0136, lwd=2, lty=2, 
        col="darkred", legend="Genome-wide proportion\nof coding sequence", 
        cex=1.4,
        bty="n", xpd=NA)

a <- tapply(1:nrow(dat), dat$nrohNorth, function(x) sum(dat$length[x] * dat$cds_proportion[x]) / sum(dat$length[x]))

bppos <- barplot(a, names.arg=names(a), ylab="",
                 main= "", cex.lab=1.5, cex.main=1.5, cex.names=1.2, xlab="", xaxt="n", yaxt="n",
                 ylim=ylims,
                 col=saola_colors2["Northern"])
axis(3, at=bppos, labels=names(a), tick=F, line=-1, cex.axis=cexaxis)
abline(h=sum(as.numeric(dat$length * dat$cds_proportion)) / sum(as.numeric(dat$length)), lwd=2, lty=2, col="darkred")

a <- tapply(1:nrow(dat), dat$nrohSouth, function(x) sum(dat$length[x] * dat$cds_proportion[x]) / sum(dat$length[x]))

bppos <- barplot(a, names.arg=names(a), ylab="",
                 main= "", cex.lab=1.5, cex.main=1.5, cex.names=1.2, xlab="", xaxt="n", yaxt="n",
                 ylim=ylims,
                 col=saola_colors2["Southern"])
axis(3, at=bppos, labels=names(a), tick=F, line=-1, cex.axis=cexaxis)
abline(h=sum(as.numeric(dat$length * dat$cds_proportion)) / sum(as.numeric(dat$length)), lwd=2, lty=2, col="darkred")
dev.off()


