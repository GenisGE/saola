f <- "/home/lpq293/mypopgen/saola/analyses/stairway_plot/results/saolaNorthern.mu4e-9.folded/saola_north.final.summary"
res <- read.table(f, h=T)

mu <- 4e-09
g <- 6


ylims <- c(0, max(res$Ne_97.5.))
xlims <- c(10, 1.5e5)

outpng <- "/home/lpq293/mypopgen/saola/paperplots2/figure3/saola_north_stairwayplot_folded.png"

saola_col <- "#86ABCB"

bitmap(outpng, h=3, w=8, res=300)

# MAKE SO Y AXIS CONSISTENT WITH PSMC
ymax <- 2e4
ylim <- c(0,ymax)
par(mar=c(5,5,4,4)+0.1, oma=c(0,1,0,0))

plot(x=res$year, y=res$Ne_median, type="l",
     xlab="",
     ylab="",
     ylim = ylim, xlim=xlims, cex.lab=2, lwd=5,
     log="x", xaxt="n", col=saola_col,cex.axis=2.5)

title(xlab=sprintf("Years ago (mu=%.2e, g=%.1f)",mu, g), cex.lab=3, line=3.5, xpd=NA)
title(ylab="Effective population size", cex.lab=3, line=3.5, xpd=NA)

lgm <- c(1.9e4, 2.65e4)

rect(xleft=lgm[1], xright = lgm[2], ybottom=-1e3, ytop=ylim[2] * 1.5,
     col="lightgrey", border=NA, xpd=FALSE)
text(y=ymax*0.9, x=2.25e4, labels="LGM", cex=3, font=2)

lines(x=res$year, y=res$Ne_median, lwd=5, col=saola_col)
lines(x=res$year, y=res$Ne_2.5., lty=2, lwd=3, col=saola_col)
lines(x=res$year, y=res$Ne_97.5., lty=2, lwd=3, col=saola_col)


xpos <- c(10, seq(100, 1000, by=500), seq(1000, 1e4, by=5e3))
xpos <- c(seq(1e1,1e2,by=1e1),
          seq(1e2,1e3,by=1e2),
          seq(1e3,1e4,by=1e3),
          seq(1e4,1e5,by=1e4),
          1e5)
xlabs <- c(5e1, 1e2, 5e2, 1e3,5e3,1e4,5e4,1e5)
xlabs <- c(1e1,2e1,5e1, 1e2,2e2,5e2, 1e3,2e3,5e3, 1e4,2e4,5e4, 1e5)
axis(1, at=xpos, labels=F)
axis(1, at=xlabs, tick=F, labels = paste(xlabs/1e3, "kya"), cex.axis=2.5)

#axis(1, at=xpos, tick=F, labels = xpos/1000, cex.axis=1.2)
dev.off()



outpdf <- "/home/lpq293/mypopgen/saola/paperplots2/figure3/saola_north_stairwayplot_folded.pdf"
pdf(outpdf, h=6, w=16)

# MAKE SO Y AXIS CONSISTENT WITH PSMC
ymax <- 2e4
ylim <- c(0,ymax)
par(mar=c(5,5,4,4)+0.1, oma=c(0,1,0,0))

plot(x=res$year, y=res$Ne_median, type="l",
     xlab="",
     ylab="",
     ylim = ylim, xlim=xlims, cex.lab=2, lwd=5,
     log="x", xaxt="n", col=saola_col,cex.axis=2.5)

title(xlab=sprintf("Years ago (mu=%.2e, g=%.1f)",mu, g), cex.lab=3, line=3.5, xpd=NA)
title(ylab="Effective population size", cex.lab=3, line=3.5, xpd=NA)

lgm <- c(1.9e4, 2.65e4)

rect(xleft=lgm[1], xright = lgm[2], ybottom=-1e3, ytop=ylim[2] * 1.5,
     col="lightgrey", border=NA, xpd=FALSE)
text(y=ymax*0.9, x=2.25e4, labels="LGM", cex=3, font=2)

lines(x=res$year, y=res$Ne_median, lwd=5, col=saola_col)
lines(x=res$year, y=res$Ne_2.5., lty=2, lwd=3, col=saola_col)
lines(x=res$year, y=res$Ne_97.5., lty=2, lwd=3, col=saola_col)


xpos <- c(10, seq(100, 1000, by=500), seq(1000, 1e4, by=5e3))
xpos <- c(seq(1e1,1e2,by=1e1),
          seq(1e2,1e3,by=1e2),
          seq(1e3,1e4,by=1e3),
          seq(1e4,1e5,by=1e4),
          1e5)
xlabs <- c(5e1, 1e2, 5e2, 1e3,5e3,1e4,5e4,1e5)
xlabs <- c(1e1,2e1,5e1, 1e2,2e2,5e2, 1e3,2e3,5e3, 1e4,2e4,5e4, 1e5)
axis(1, at=xpos, labels=F)
axis(1, at=xlabs, tick=F, labels = paste(xlabs/1e3, "kya"), cex.axis=2.5)

#axis(1, at=xpos, tick=F, labels = xpos/1000, cex.axis=1.2)
dev.off()
