
covfile <- "/home/lpq293/mypopgen/saola/analyses/beagle_pcangsd/results_finalnorelated/pcangsd/out.cov"
popfile <- "/home/lpq293/mypopgen/saola/analyses/beagle_pcangsd/results_finalnorelated/info/pop.list"
bamlist <- "/home/lpq293/mypopgen/saola/analyses/beagle_pcangsd/results_finalnorelated/info/bams.list"

covfile <- "/home/lpq293/mypopgen/saola/analyses/beagle_pcangsd/results_mergedseq/pcangsd/out.cov"
popfile <- "/home/lpq293/mypopgen/saola/analyses/beagle_pcangsd/results_mergedseq/info/pop.list"
bamlist <- "/home/lpq293/mypopgen/saola/analyses/beagle_pcangsd/results_mergedseq/info/bams.list"

cov <- as.matrix(read.table(covfile))
ei <- eigen(cov)
e <- ei$vectors
v <- ei$values
vars <- v/sum(v)*100

pop <- scan(popfile, what="dfas")
pop <- gsub("Central", "Southern", pop)
npop <- length(unique(pop))

bams <- scan(bamlist, what="d")
ids <- gsub(".saolaRefN.bam", "", basename(bams))

saola_colors <- c("Northern" = "#86ABCB",
                  "Southern" = "#FFB300")



outpng <- "/home/lpq293/mypopgen/saola/paperplots2/figure2/pca.png"
bitmap(outpng, width=2.5, height=2.5, res=300)
par(mar=c(4,5,2,1.5))
#xlegend <- max(e[,1]) * 1.1
#ylegend <- max(e[,2])         
plot(e[,1], e[,2], pch=21,
         cex=3.5, bg=saola_colors[pop], 
         xlab=paste0("PC 1 (",round(vars[1], 2),"%)"), 
         ylab=paste0("PC 2 (",round(vars[2], 2),"%)"),
         cex.lab=2.5, cex.axis=1.9, xpd=NA)
#legend("topright",
#       legend=unique(pop), pt.bg=saola_colors,
#       pch=21, pt.cex=2,bty='n', xpd=NA, cex=2)
dev.off()



outpdf <- "/home/lpq293/mypopgen/saola/paperplots2/figure2/pca.pdf"
pdf(outpdf, width=4.4, height=4.4)#, res=300)
par(mar=c(4,5,2,1.5))
#xlegend <- max(e[,1]) * 1.1
#ylegend <- max(e[,2])         
plot(e[,1], e[,2], pch=21,
         cex=3.5, bg=saola_colors[pop], 
         xlab=paste0("PC 1 (",round(vars[1], 2),"%)"), 
         ylab=paste0("PC 2 (",round(vars[2], 2),"%)"),
         cex.lab=2.5, cex.axis=1.9, xpd=NA)
#legend("topright",
#       legend=unique(pop), pt.bg=saola_colors,
#       pch=21, pt.cex=2,bty='n', xpd=NA, cex=2)
dev.off()

