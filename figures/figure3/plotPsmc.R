source("/home/lpq293/mypopgen/saola/info/loadPopInfo.R")

continuousPal <- function(x, cols=c("white", "darkblue"), n=100, maxz=20){

    colpal <- colorRampPalette(cols)
    wholepal <- colpal(n)

    x[x>maxz] <- maxz

    mapvals <- seq(range(x)[1], range(x)[2], length.out=n)
    return(wholepal[sapply(x, function(x) which.min(abs(x-mapvals)))])
}



names(pop) <- id
pop <- gsub("Central", "Southern", pop)
saola_colors <- c("Northern" = "#86ABCB",
                  "Southern" = "#FFB300")


mu <- 1.2e-08
g <- 6

source("readPSMC.R")


outpng <- "saolaPsmcCrossCoal.png"
bitmap(outpng, w=8, h=4, res=300)
ymax <- 2e4 # ARBITRARY THRESHOLD BECAUSE PLOTS USUALY LOOK NICE WITH THIS
#par(mar=c(5,5,4,4)+0.1, oma=c(0,0,4,0)) # HORIZONTAL TOP LEGEND
par(mar=c(5,5,4,4)+0.1, oma=c(0,0.5,2,4)) # VERTICAL RIGHT LEGEND
plot(type='l', x=res_psmc[[1]]$YearsAgo, log='x', y=res_psmc[[1]]$Ne,
     col = saola_colors[pop[names(res_psmc[1])]], lwd=5,
     xlab="",
     ylab="",cex.lab=2, xaxt='n',
     xlim=c(2 * 10^3, 1.5 * 10^6), ylim=c(0, ymax), cex.axis=2.5,
     main="")

title(xlab=sprintf("Years ago (mu=%.2e, g=%.1f)",mu, g), cex.lab=3, line=3.5, xpd=NA)
title(ylab="Effective population size", cex.lab=3, line=3.5, xpd=NA)


lgm <- c(1.9e4, 2.65e4)

rect(xleft=lgm[1], xright = lgm[2], ybottom=-1e3, ytop=ymax,
     col="lightgrey", border=NA)
text(y=ymax*0.7, x=2.25e4, labels="LGM", cex=3, font=2)

#xpos <- seq(1.5 * 10^4, 1.5 * 10^6,by=10000)
xpos <- c(seq(1e3, 1e4, by=1e3), seq(1e4,1e5,by=1e4), seq(1e5, 1e6, by=1e5))
xlabs <- c(2e3, 1e4, 1e5, 5e5,1e6)
axis(1, at=xpos, labels=F)
axis(1, at=xlabs, tick=F, labels= paste(as.integer(xlabs/1e3), "kya"), cex.axis=2.5)
for(i in 1:length(res_psmc)) lines(x=res_psmc[[i]]$YearsAgo, y=res_psmc[[i]]$Ne,
                col = saola_colors[pop[names(res_psmc[i])]],
               lwd=5)

legend(legend=names(saola_colors),
       col=saola_colors,
       lwd=5,border=NA, bty="n",cex=3.3, xpd=NA,
       y=8e3, x=1.5e5)

# plot crosscoalescence
crosscoalpalcols <- c("white", "darkred")


rect(xleft=crosscoaltimes[-length(crosscoaltimes)],
     xright=crosscoaltimes[-1],
     ybottom = ymax - ymax * 0.15,
     ytop = ymax * 2,
     border=NA, col=continuousPal(crosscoalrate, crosscoalpalcols),
     xpd=F)

segments(x0=1e3, x1=2e6, y0=ymax - ymax * 0.15)

## add legend
n <- 100
colpal <- colorRampPalette(crosscoalpalcols)
wholepal <- colpal(n)

if(TRUE){ # THIS MAKES COLORSCALE LEGEND VERTICAL TO THE RIGHT
    
    rasterImage(rev(wholepal),
                xleft=2.15e6, xright=3e6,
                ybottom=1.2e4, ytop=2.2e4,
                xpd=NA)
    rect(xleft=2.15e6, xright=3e6,
         ybottom=1.2e4, ytop=2.2e4, xpd=NA)
    #text(x=2e5, y=2.65e4, labels="Relative cross-coalescence rate", xpd=NA, cex = 1.8, font=2)
    text(x= 3.7e6, y=c(1.2e4, 2.2e4), labels=signif(range(crosscoalrate),1), xpd=NA, cex = 3)


}


if(FALSE){ # THIS MAKES COLORSCALE LEGEND HORIZONTAL IN TOP
    rasterImage(t(wholepal),
                xleft=1e5, xright=5e5,
                ybottom=2.25e4, ytop=2.55e4,
                xpd=NA)
    rect(xleft=1e5, xright=5e5,
         ybottom=2.25e4, ytop=2.55e4, xpd=NA)
    #text(x=2e5, y=2.65e4, labels="Relative cross-coalescence rate", xpd=NA, cex = 1.8, font=2)
    text(x=c(1e5, 5e5), y=2.17e4, labels=signif(range(crosscoalrate),1), xpd=NA, cex = 2)
}

#text(x=0.8, y = c(0.25,0.5, 0.75),
#     labels = c(-max(abs(min_z),abs(max_z)), 0, max(abs(min_z),abs(max_z))),
#     cex=cex.legend,xpd=NA)

# this just redraws the box so it is not removed by crosscoalescence rect
par(new=TRUE)
plot(1, log="x", type="n",
     xlab="", yaxt="n",
     ylab="", xaxt="n",
     xlim=c(2 * 10^3, 1.5 * 10^6), ylim=c(0, ymax),
     main="")

text(x=2e3, y=ymax + ymax * 0.1,
     labels="Relative cross coalescence rate:",
     col="black", cex=2.5,font=2, xpd=NA, adj=0)
text(x=1.9e3, y=ymax - ymax * 0.27,
     labels="Effective\npopulation size:",
     col="black", cex=2.5,font=2, xpd=NA, adj=0)

dev.off()




outpdf <- "saolaPsmcCrossCoal.pdf"
pdf(outpdf, w=16, h=8)#, res=300)
ymax <- 2e4 # ARBITRARY THRESHOLD BECAUSE PLOTS USUALY LOOK NICE WITH THIS
#par(mar=c(5,5,4,4)+0.1, oma=c(0,0,4,0)) # HORIZONTAL TOP LEGEND
par(mar=c(5,5,4,4)+0.1, oma=c(0,0.5,2,4)) # VERTICAL RIGHT LEGEND
plot(type='l', x=res_psmc[[1]]$YearsAgo, log='x', y=res_psmc[[1]]$Ne,
     col = saola_colors[pop[names(res_psmc[1])]], lwd=5,
     xlab="",
     ylab="",cex.lab=2, xaxt='n',
     xlim=c(2 * 10^3, 1.5 * 10^6), ylim=c(0, ymax), cex.axis=2.5,
     main="")

title(xlab=sprintf("Years ago (mu=%.2e, g=%.1f)",mu, g), cex.lab=3, line=3.5, xpd=NA)
title(ylab="Effective population size", cex.lab=3, line=3.5, xpd=NA)


lgm <- c(1.9e4, 2.65e4)

rect(xleft=lgm[1], xright = lgm[2], ybottom=-1e3, ytop=ymax,
     col="lightgrey", border=NA)
text(y=ymax*0.7, x=2.25e4, labels="LGM", cex=3, font=2)

#xpos <- seq(1.5 * 10^4, 1.5 * 10^6,by=10000)
xpos <- c(seq(1e3, 1e4, by=1e3), seq(1e4,1e5,by=1e4), seq(1e5, 1e6, by=1e5))
xlabs <- c(2e3, 1e4, 1e5, 5e5,1e6)
axis(1, at=xpos, labels=F)
axis(1, at=xlabs, tick=F, labels= paste(as.integer(xlabs/1e3), "kya"), cex.axis=2.5)
for(i in 1:length(res_psmc)) lines(x=res_psmc[[i]]$YearsAgo, y=res_psmc[[i]]$Ne,  
                    col = saola_colors[pop[names(res_psmc[i])]],
                    lwd=5)

legend(legend=names(saola_colors),
       col=saola_colors,
       lwd=5,border=NA, bty="n",cex=3.3, xpd=NA,
       y=8e3, x=1.5e5)

# plot crosscoalescence
crosscoalpalcols <- c("white", "darkred")


rect(xleft=crosscoaltimes[-length(crosscoaltimes)],
     xright=crosscoaltimes[-1],
     ybottom = ymax - ymax * 0.15,
     ytop = ymax * 2,
     border=NA, col=continuousPal(crosscoalrate, crosscoalpalcols),
     xpd=F)

segments(x0=1e3, x1=2e6, y0=ymax - ymax * 0.15)

## add legend
n <- 100
colpal <- colorRampPalette(crosscoalpalcols)
wholepal <- colpal(n)

if(TRUE){ # THIS MAKES COLORSCALE LEGEND VERTICAL TO THE RIGHT
    
    rasterImage(rev(wholepal),
                xleft=2.15e6, xright=3e6,
                ybottom=1.2e4, ytop=2.2e4,
                xpd=NA)
    rect(xleft=2.15e6, xright=3e6,
         ybottom=1.2e4, ytop=2.2e4, xpd=NA)
    #text(x=2e5, y=2.65e4, labels="Relative cross-coalescence rate", xpd=NA, cex = 1.8, font=2)
    text(x= 3.7e6, y=c(1.2e4, 2.2e4), labels=signif(range(crosscoalrate),1), xpd=NA, cex = 3)


}


if(FALSE){ # THIS MAKES COLORSCALE LEGEND HORIZONTAL IN TOP
    rasterImage(t(wholepal),
                xleft=1e5, xright=5e5,
                ybottom=2.25e4, ytop=2.55e4,
                xpd=NA)
    rect(xleft=1e5, xright=5e5,
         ybottom=2.25e4, ytop=2.55e4, xpd=NA)
    #text(x=2e5, y=2.65e4, labels="Relative cross-coalescence rate", xpd=NA, cex = 1.8, font=2)
    text(x=c(1e5, 5e5), y=2.17e4, labels=signif(range(crosscoalrate),1), xpd=NA, cex = 2)
}

#text(x=0.8, y = c(0.25,0.5, 0.75),
#     labels = c(-max(abs(min_z),abs(max_z)), 0, max(abs(min_z),abs(max_z))),
#     cex=cex.legend,xpd=NA)

# this just redraws the box so it is not removed by crosscoalescence rect
par(new=TRUE)
plot(1, log="x", type="n",
     xlab="", yaxt="n",
     ylab="", xaxt="n",
     xlim=c(2 * 10^3, 1.5 * 10^6), ylim=c(0, ymax),
     main="")

text(x=2e3, y=ymax + ymax * 0.1,
     labels="Relative cross coalescence rate:",
     col="black", cex=2.5,font=2, xpd=NA, adj=0)
text(x=1.9e3, y=ymax - ymax * 0.27,
     labels="Effective\npopulation size:",
     col="black", cex=2.5,font=2, xpd=NA, adj=0)


dev.off()
