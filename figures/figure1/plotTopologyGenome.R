

if(FALSE){
getwins <- function(d, w=5e6, chrord=chrord){

    maxes <- tapply(d$pos, d$chr, max)[chrord]
    winbreaks <- list()
    
    chrs <- unique(d$chr)

    chrss <- c()
    midss <- c()
    topo0 <- c()
    topo1 <- c()
    topo2 <- c()
    topo3 <- c()
   
    
    for(c in chrs){

        winbreaks <- seq(0, maxes[c]+w, by=w)
        wins <- cut(d$pos[d$chr==c], winbreaks)

        chrss <- c(chrss, rep(c, length(winbreaks) - 1))
        midss <- c(midss, winbreaks[-1] - w/2)
        topo0 <- c(topo0, tapply(d$topos[d$chr==c], wins,  function(x) mean(x == 0)))
        topo1 <- c(topo1, tapply(d$topos[d$chr==c], wins,  function(x) mean(x == 1)))
        topo2 <- c(topo2, tapply(d$topos[d$chr==c], wins,  function(x) mean(x == 2)))
        topo3 <- c(topo3, tapply(d$topos[d$chr==c], wins,  function(x) mean(x == 3)))
        
    }

    dwin <- data.frame(chr=chrss, midpos = midss,
                       topo0 = topo0, topo1=topo1, topo2=topo2, topo3=topo3,
                       stringsAsFactors=F)
    dwin <- dwin[!is.na(dwin$topo0),]
    return(dwin)

}
}

ordchr <- function(x, c){
    idx <- which(x$chr==c)
    ord <- order(x$pos[idx])
    idx[ord]
}



slidingwins <- function(d, window.size=1e6, step.size=1e5){
    
        
    chroms <- unique(d$chr)
    mid.poss <- c()
    win.chroms <- c()
    
    for(c in chroms){

        maxchrom <- max(d$pos[d$chr==c])
        starts <- seq(0, maxchrom, by=step.size)
        mid.poss <- c(mid.poss, starts +  window.size / 2)
        
        win.chroms <- c(win.chroms, rep(c, length(starts)))
    }

    n <- length(mid.poss)
    nwins <- integer(n)
    topo0 <- numeric(n)
    topo1 <- numeric(n)
    topo2 <- numeric(n)
    topo3 <- numeric(n)

    # estiamted proproitn of topologies per window
    for(i in 1:n){

        c <- win.chroms[i]
        start <- mid.poss[i] - window.size / 2
        end <- mid.poss[i] + window.size / 2
        k <- d$chr == c & d$pos > start & d$pos < end
        nwins[i] <- sum(k)
        topo0[i] <- mean(d$topo[k] == 0)
        topo1[i] <- mean(d$topo[k] == 1)
        topo2[i] <- mean(d$topo[k] == 2)
        topo3[i] <- mean(d$topo[k] == 3)


    }

    outd <- data.frame(chr=win.chroms, pos = mid.poss, nwins = nwins,
                       topo0 = topo0, topo1=topo1, topo2=topo2, topo3=topo3,
                       stringsAsFactors=F)
    return(outd)
}

f <- "files/topo.all.gz"

d <- read.table(f, stringsAsFactors=F, h=F)

chrpos <- do.call("rbind", strsplit(d$V1, "_"))
topos <- data.frame(chr=chrpos[,1], pos=as.integer(chrpos[,2]), topos = d$V2, stringsAsFactors=F)

chrord <- c(paste0("chr", 1:29), "chrX")

ord <- unlist(sapply(chrord, ordchr, x=topos))
topos <- topos[ord,]

if(FALSE){
    ## get some statistics of topology for tree
    ## genomewide percentage of topologies
    table(topos$topos)[1]/sum(table(topos$topos)) 
    table(topos$topos)[2]/sum(table(topos$topos))
    table(topos$topos)[3]/sum(table(topos$topos))

    xchrom <- topos$chr == "chrX"
    ## autosomal percentage
    table(topos$topos[!xchrom])[1]/sum(table(topos$topos[!xchrom]))
    table(topos$topos[!xchrom])[2]/sum(table(topos$topos[!xchrom]))
    table(topos$topos[!xchrom])[3]/sum(table(topos$topos[!xchrom]))

    ## x chromosome percentage
    table(topos$topos[xchrom])[1]/sum(table(topos$topos[xchrom]))
    table(topos$topos[xchrom])[2]/sum(table(topos$topos[xchrom]))
    table(topos$topos[xchrom])[3]/sum(table(topos$topos[xchrom]))
}

wins1 <- slidingwins(topos)
wins1 <- wins1[!wins1$nwins<median(wins1$nwins),]

#barplot(t(wins1[,4:7]), border=NA, space=0, xaxt="n")

wins2 <- slidingwins(topos,  window.size=5e6, step.size=5e5)
k <- !wins2$nwins<median(wins2$nwins) | (wins2$chr == "chrX" & wins2$nwins > median(wins2$nwins[wins2$chr == "chrX"]))

wins2 <- wins2[k,]


outpng <- "/home/lpq293/mypopgen/saola/paperplots2/figure1/topologiesFrac.png"

bitmap(outpng, h=3,w=12,res=300)
colpal <- c("mediumpurple3", "indianred", "dodgerblue2")

par(oma=c(4,5,0,13))

barplot(t(wins2[,4:7]), border=NA, space=0, xaxt="n", col=colpal,
        ylab="", cex.lab=2.5, cex.axis=2.4, xpd=NA)
title(ylab="Fraction of blocks\nwith topology", cex.lab=2.5, line=3.5, xpd=NA)

chrbreaks <- c(0,cumsum(sapply(unique(wins2$chr), function(x){sum(wins2$chr==x)})))
abline(v=chrbreaks, col="white", lwd=1)
segments(x0=chrbreaks[-length(chrbreaks)], x1=chrbreaks[-1], y0=-c(0.14, 0.2), xpd=NA, lwd=5)

text(x = sort(tapply(1:nrow(wins2),wins2$chr,mean)),
     y = -c(0.07, 0.28),
     labels=gsub("chr", "", unique(wins2$chr)),
     xpd=NA, cex=2.5)
text(x=nrow(wins2)/2, y=-0.44, labels="Chromosome in cattle", xpd=NA, cex=3, font=2)



xat <- c(nrow(wins2) * 1.02,
            nrow(wins2) * 1.05,
            nrow(wins2) * 1.08)

yat <- c(0.15, 0.075, 0, 0.125, 0.25)

ymove <- c(-0.3, 0.2, 0.7)

segments(x0=xat[1], x1=xat[3],
            y0=yat[1] + ymove, y1=yat[3] + ymove,
             xpd=NA, lwd=4,
             col=colpal)
segments(x0=xat[1], x1=xat[3],
            y0=yat[1] + ymove, y1=yat[5] + ymove,
             xpd=NA, lwd=4,
             col=colpal)
segments(x0=xat[2], x1=xat[3],
            y0=yat[2] + ymove, y1=yat[4] + ymove,
             xpd=NA, lwd=4,
             col=colpal)

xtext <- nrow(wins2) * 1.085
ytext <- c(0.25, 0.125, 0)
topo1 <- c("Saola", "Water Buffalo", "Cattle")
topo2 <- c("Cattle", "Saola", "Water buffalo") 
topo3 <- c("Water buffalo", "Saola", "Cattle")

text(x=nrow(wins2) * 1.1, y=1.1, labels="Topologies",
        cex=2.5, font=2, xpd=NA)
text(x=xtext, y=ytext + ymove[1], labels=topo1, 
    col=colpal[1], 
    font=2,cex=2, xpd=NA, adj=0)
text(x=xtext, y=ytext + ymove[2], labels=topo2, 
    col=colpal[2], 
    font=2,cex=2, xpd=NA, adj=0)

text(x=xtext, y=ytext + ymove[3], labels=topo3, 
    col=colpal[3], 
    font=2,cex=2, xpd=NA, adj=0)

#text(y=0.9, x= nrow(wins2) * 1.06, labels="Topologies", cex=1.7, xpd=NA, font=2)
#text(y=c(0.65, 0.5, 0.35)[3:1] + 0.1,
#     x= nrow(wins2) * 1.1,
#     labels=toposlab, col=colpal, xpd=NA, font=2, cex=1.5)
dev.off()




outpdf <- "/home/lpq293/mypopgen/saola/paperplots2/figure1/topologiesFrac.pdf"
pdf(outpdf, h=6,w=24)
colpal <- c("mediumpurple3", "indianred", "dodgerblue2")

par(oma=c(4,5,0,13))

barplot(t(wins2[,4:7]), border=NA, space=0, xaxt="n", col=colpal,
        ylab="", cex.lab=2.5, cex.axis=2.4, xpd=NA)
title(ylab="Fraction of blocks\nwith topology", cex.lab=2.5, line=3.5, xpd=NA)

chrbreaks <- c(0,cumsum(sapply(unique(wins2$chr), function(x){sum(wins2$chr==x)})))
abline(v=chrbreaks, col="white", lwd=1)
segments(x0=chrbreaks[-length(chrbreaks)], x1=chrbreaks[-1], y0=-c(0.14, 0.2), xpd=NA, lwd=5)

text(x = sort(tapply(1:nrow(wins2),wins2$chr,mean)),
     y = -c(0.07, 0.28),
     labels=gsub("chr", "", unique(wins2$chr)),
     xpd=NA, cex=2.5)
text(x=nrow(wins2)/2, y=-0.44, labels="Chromosome in cattle", xpd=NA, cex=3, font=2)



xat <- c(nrow(wins2) * 1.02,
            nrow(wins2) * 1.05,
            nrow(wins2) * 1.08)

yat <- c(0.15, 0.075, 0, 0.125, 0.25)

ymove <- c(-0.3, 0.2, 0.7)

segments(x0=xat[1], x1=xat[3],
            y0=yat[1] + ymove, y1=yat[3] + ymove,
             xpd=NA, lwd=4,
             col=colpal)
segments(x0=xat[1], x1=xat[3],
            y0=yat[1] + ymove, y1=yat[5] + ymove,
             xpd=NA, lwd=4,
             col=colpal)
segments(x0=xat[2], x1=xat[3],
            y0=yat[2] + ymove, y1=yat[4] + ymove,
             xpd=NA, lwd=4,
             col=colpal)

xtext <- nrow(wins2) * 1.085
ytext <- c(0.25, 0.125, 0)
topo1 <- c("Saola", "Water Buffalo", "Cattle")
topo2 <- c("Cattle", "Saola", "Water buffalo") 
topo3 <- c("Water buffalo", "Saola", "Cattle")

text(x=nrow(wins2) * 1.1, y=1.1, labels="Topologies",
        cex=2.5, font=2, xpd=NA)
text(x=xtext, y=ytext + ymove[1], labels=topo1, 
    col=colpal[1], 
    font=2,cex=2, xpd=NA, adj=0)
text(x=xtext, y=ytext + ymove[2], labels=topo2, 
    col=colpal[2], 
    font=2,cex=2, xpd=NA, adj=0)

text(x=xtext, y=ytext + ymove[3], labels=topo3, 
    col=colpal[3], 
    font=2,cex=2, xpd=NA, adj=0)

#text(y=0.9, x= nrow(wins2) * 1.06, labels="Topologies", cex=1.7, xpd=NA, font=2)
#text(y=c(0.65, 0.5, 0.35)[3:1] + 0.1,
#     x= nrow(wins2) * 1.1,
#     labels=toposlab, col=colpal, xpd=NA, font=2, cex=1.5)
dev.off()

