source("/home/lpq293/mypopgen/saola/info/loadPopInfo.R")
source("/home/lpq293/github/evalAdmix/visFuns.R")
source("/home/lpq293/mypopgen/saola/paperplots/figure2/aaplot_gge2.R")
#source("/home/albrecht/github/evalAdmix/aaplot.R")

qfile <- "/home/lpq293/mypopgen/saola/analyses/ngsadmix/results_all26/2/admixNotransSaolaall26.2.3.qopt_conv"


ids <- df$paperIDs[is.na(df$exclude_reason)]
newids <- df$newPaperIDs[is.na(df$exclude_reason)]
# make N go first than C
num <- as.integer(gsub("[NC]", "", ids))
p <- gsub("[0-9]", "", ids)
ord <- order(ifelse(p=="N",1,2), num)

## read in admixture proporitons and correlation of residuals
q <- as.matrix(read.table(qfile))

## admix and evaladmix have different sample ord. make consistent than depth
admix_ids <- gsub(".saolaRefN.bam", "", basename(read.table("/home/lpq293/mypopgen/saola/analyses/beagle_pcangsd/saola_26_samples_mergedseq_generalpop_pop.info",h=F, stringsAsFactors=F)$V2))
ord2 <- sapply(1:length(admix_ids), function(x) which(admix_ids == id[ord][x]))




### read in depths
indir <- "/home/lpq293/mypopgen/saola/analyses/coverage/results/samples"
depths <- list()

for(f in list.files(indir, "depth", full.names=T)[gsub(".depth", "",list.files(indir, "depth")) %in% id]) depths[[gsub(".depth","", basename(f))]] <- scan(f, what=2)[3]

d <- unlist(depths)[id]

outpng <- "/home/lpq293/mypopgen/saola/paperplots2/figure2/depthAdmix.png"
#bitmap(outpng, h=5,w=5,res=300)
png(outpng, h=6,w=10, units="in",res=300)

par(mfrow=c(2,1), mar=c(1.5,4,1,1),oma=c(0,1,0,0))
x <- barplot(log10(d + 1)[ord], xlab="",  yaxt="n",
        ylab="Sequencing depth",
        names.arg="",#ids[ord],
        col=saola_colors[pop[ord]],
        cex.lab=1.5, las=2, xpd=NA)

axis(1, at=x[,1], labels=newids[ord], cex.axis=1.5, las=2, tick=F, line=-0.75)
xpos <- c(0,1,2.5, 5,10,20,40)
axis(2, at=log10(xpos+1), labels=xpos,xpd=NA)
abline(h=log10(xpos + 1)[-1], lty=2, lwd=0.3)
barplot(log10(d + 1)[ord], xlab="",  yaxt="n",
        ylab="",
        names.arg="",#ids[ord],
        col=saola_colors[pop[ord]], add=T)

par(mar=c(1,4,3,1))
barplot(t(q)[,ord2], col=saola_colors[2:1], space=0, border=NA, cex.axis=1.2,cex.lab=1.5,
        ylab="Admixture proportions", xlab="",xpd=NA)
abline(v=1:nrow(q), col="white", lwd=0.2)

abline(v=cumsum(sapply(unique(pop[ord])[1],function(x){sum(pop[ord]==x)})),col=1,lwd=1.2)

dev.off()



outpdf <- "/home/lpq293/mypopgen/saola/paperplots2/figure2/depthAdmix.pdf"
#bitmap(outpng, h=5,w=5,res=300)
pdf(outpdf, h=6,w=10)#, units="in")#,res=300)

par(mfrow=c(2,1), mar=c(1.5,4,1,1),oma=c(0,1,0,0))
x <- barplot(log10(d + 1)[ord], xlab="",  yaxt="n",
        ylab="Sequencing depth",
        names.arg="",#ids[ord],
        col=saola_colors[pop[ord]],
        cex.lab=1.5, las=2, xpd=NA)

axis(1, at=x[,1], labels=newids[ord], cex.axis=1.5, las=2, tick=F, line=-0.75)
xpos <- c(0,1,2.5, 5,10,20,40)
axis(2, at=log10(xpos+1), labels=xpos,xpd=NA)
abline(h=log10(xpos + 1)[-1], lty=2, lwd=0.3)
barplot(log10(d + 1)[ord], xlab="",  yaxt="n",
        ylab="",
        names.arg="",#ids[ord],
        col=saola_colors[pop[ord]], add=T)

par(mar=c(1,4,3,1))
barplot(t(q)[,ord2], col=saola_colors[2:1], space=0, border=NA, cex.axis=1.2,cex.lab=1.5,
        ylab="Admixture proportions", xlab="",xpd=NA)
abline(v=1:nrow(q), col="white", lwd=0.2)

abline(v=cumsum(sapply(unique(pop[ord])[1],function(x){sum(pop[ord]==x)})),col=1,lwd=1.2)

dev.off()
