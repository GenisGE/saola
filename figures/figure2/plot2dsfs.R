# contains function "sfs_hudsons_fst" to calcualte Hudson's fst from 2dsfs
source("/home/lpq293/github/QCSeq/gt/scripts/sfs.R")


fold2dsfs <- function(sfsmat){

    outsfs <- matrix(0, nrow=nrow(sfsmat),ncol=ncol(sfsmat))
    n_half <- (nrow(sfsmat) + ncol(sfsmat) - 2) / 2
    
    for(i in 1:nrow(sfsmat)){
        for(j in 1:ncol(sfsmat)){
            if((i+j-2) < n_half){
                outsfs[i,j] <- outsfs[i,j] + sfsmat[i,j]
            } else if ((i+j) < n_half) {
                outsfs[nrow(outsfs) - i,ncol(outsfs) - j] <- outsfs[i,j] + sfsmat[i,j]
            } else if ((i+j) == n_half) {
                outsfs[i,j] <- outsfs[i,j] + sfsmat[i,j] / 2
                outsfs[nrow(outsfs) - i,ncol(outsfs) - j] <- outsfs[i,j] + sfsmat[i,j] / 2
            }
        }
    }
    return(outsfs)
}



f <- "/home/lpq293/mypopgen/saola/analyses/sfs_unfolded/results/notrans/2dsfs/Central_Northern2.sfs"

sfs <- scan(f, what=.3, skip=1)
sfsh <- scan(f, what="d", nlines=1, sep="=")
dim <- as.numeric(strsplit(gsub("[<>]","",sfsh[2]), split="/")[[1]])

sfsmat <- matrix(sfs, nrow=dim[2], ncol=dim[1])
foldedsfs <- fold2dsfs(sfsmat)
foldedsfs2 <- foldedsfs
foldedsfs2[1,1] <- 0


colpal <- colorRampPalette(colors = c("lightgrey", "darkred"), space="Lab")(50)
ncol <- 9
colpal <- RColorBrewer::brewer.pal("YlGnBu", n=ncol)

outpng <- "/home/lpq293/mypopgen/saola/paperplots2/figure2/northernSouthern2dsfs.png"
fst <- sfs_hudsons_fst(foldedsfs2)


bitmap(outpng, w=6,h=5.5, res=300)
par(oma=c(0.5,2.2,0,9))

toplot <- log10(t(foldedsfs2/sum(foldedsfs)))
toplot[is.infinite(toplot)] <- NA
freqvals <- range(toplot, na.rm=T)

image(0:(ncol(foldedsfs2)-1), 0:(nrow(foldedsfs2)-1), toplot, 
    col = colpal, axes=FALSE,
    ylab="", 
    xlab="", main="",
      cex.lab=3.5, xpd=NA)

title(ylab="Northern minor allele count", cex.lab=3.2, xpd=NA, line=4)
title(xlab="Southern minor allele count", cex.lab=3.2, xpd=NA, line=4)

axis(1, 0:(ncol(foldedsfs2)-1), colnames(foldedsfs2), tick=FALSE, xpd=NA, cex.axis=2.5)
axis(2, 0:(nrow(foldedsfs2)-1), rownames(foldedsfs2), las=2, tick=FALSE, xpd=NA, cex.axis=2.5)
text(x=0, y=0, labels=signif(foldedsfs[1] / sum(foldedsfs), 4), cex=1.5)

rect(ybottom=seq(-0.5, nrow(foldedsfs2)-0.5, length.out=ncol+1)[-(ncol + 1)],
     ytop = seq(-0.5, nrow(foldedsfs2)-0.5, length.out=ncol+1)[-1],
     xleft=ncol(foldedsfs2), xright=ncol(foldedsfs2) + 0.5,
     col=colpal, xpd=NA)

text(x=ncol(foldedsfs2) + 1.5,
     y=seq(-0.5, nrow(foldedsfs2)-0.5, length.out=ncol+1),
     labels=signif(10 ^ seq(freqvals[1], freqvals[2], length.out=ncol+1), 2),
     xpd=NA, cex=1.9)
text(y=nrow(foldedsfs2)+0.12, x=ncol(foldedsfs2) + 0.5, labels="Frequency", cex=2.5, xpd=NA)

text(y=nrow(foldedsfs2) * 0.7, x=ncol(foldedsfs2) * 0.5,
     labels=bquote(paste("F"["ST"]*" = ", .(round(fst, 2)))), cex=3.5)
         #paste0(expression("F"["ST"]" = "), round(fst, 2)))

dev.off()




outpdf <- "/home/lpq293/mypopgen/saola/paperplots2/figure2/northernSouthern2dsfs.pdf"
#fst <- sfs_hudsons_fst(foldedsfs2)

pdf(outpdf, w=9,h=8.5)
par(oma=c(0.5,2.2,0,9))

toplot <- log10(t(foldedsfs2/sum(foldedsfs)))
toplot[is.infinite(toplot)] <- NA
freqvals <- range(toplot, na.rm=T)

image(0:(ncol(foldedsfs2)-1), 0:(nrow(foldedsfs2)-1), toplot, 
    col = colpal, axes=FALSE,
    ylab="", 
    xlab="", main="",
      cex.lab=3.5, xpd=NA)

title(ylab="Northern minor allele count", cex.lab=3.2, xpd=NA, line=4)
title(xlab="Southern minor allele count", cex.lab=3.2, xpd=NA, line=4)

axis(1, 0:(ncol(foldedsfs2)-1), colnames(foldedsfs2), tick=FALSE, xpd=NA, cex.axis=2.5)
axis(2, 0:(nrow(foldedsfs2)-1), rownames(foldedsfs2), las=2, tick=FALSE, xpd=NA, cex.axis=2.5)
text(x=0, y=0, labels=signif(foldedsfs[1] / sum(foldedsfs), 4), cex=1.5)

rect(ybottom=seq(-0.5, nrow(foldedsfs2)-0.5, length.out=ncol+1)[-(ncol + 1)],
     ytop = seq(-0.5, nrow(foldedsfs2)-0.5, length.out=ncol+1)[-1],
     xleft=ncol(foldedsfs2), xright=ncol(foldedsfs2) + 0.5,
     col=colpal, xpd=NA)

text(x=ncol(foldedsfs2) + 1.5,
     y=seq(-0.5, nrow(foldedsfs2)-0.5, length.out=ncol+1),
     labels=signif(10 ^ seq(freqvals[1], freqvals[2], length.out=ncol+1), 2),
     xpd=NA, cex=1.9)
text(y=nrow(foldedsfs2)+0.12, x=ncol(foldedsfs2) + 0.5, labels="Frequency", cex=2.5, xpd=NA)

text(y=nrow(foldedsfs2) * 0.7, x=ncol(foldedsfs2) * 0.5,
     labels=bquote(paste("F"["ST"]*" = ", .(round(fst, 2)))), cex=3.5)
         #paste0(expression("F"["ST"]" = "), round(fst, 2)))


dev.off()
