

cexaxis <- 1.8
cexlab <- 2.5
bitmap("figure6.png", h=6,w=12, res=300)
layout(matrix(1:6, nrow=3))
source("plotPast.R")
cat("Finished first plot\n")
source("plotFuture.R")
text(y=ylims[2] * 3.47, x=-c(xlims[2]*1.34, xlims[2] * 0.13), labels=c("A", "B"), xpd=NA, cex=4, font=2)
dev.off()


cexaxis <- 1.8
cexlab <- 2.5
pdf("figure6.pdf", h=8,w=16)
layout(matrix(1:6, nrow=3))
source("plotPast.R")
cat("Finished first plot\n")
source("plotFuture.R")
text(y=ylims[2] * 3.47, x=-c(xlims[2]*1.34, xlims[2] * 0.13), labels=c("A", "B"), xpd=NA, cex=4, font=2)
dev.off()
