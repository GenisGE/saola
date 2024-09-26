source("/home/lpq293/mypopgen/saola/info/loadPopInfo.R")

f <- "/home/lpq293/mypopgen/saola/analyses/roh/plink/results/roh/roh_min1000kb_den100_win100_het1.hom" # plink .hom file with rohs
g <- "/home/lpq293/mypopgen/saola/genome_masks/v2/saola_ref/all/beds/saolaRefN_autosomal_10MBscaff.bed" # bed file with size of autosomes
#fam <- "samples_plot.info" # paht to fam file, assumed that individual ids are first column
#minlen <- as.numeric(args[5]) / 1000

autosome <- read.table(g)
autosome <- autosome[autosome$V3 > 5e6,]


autosome_len <- sum(as.numeric(autosome$V3))

roh <- read.table(f, h=T)
roh <- roh[roh$CHR %in% autosome$V1,]
#roh <- roh[roh$KB > minlen,]

getROHproportion <- function(roh, s, autosome_len){

    x <- roh[roh$IID==s,]
    d <- c(
        sum(x$KB[x$KB >= 1e3 &x$KB < 2.5e3] * 1e3),
        sum(x$KB[x$KB >= 2.5e3 & x$KB < 5e3] * 1e3),
        sum(x$KB[x$KB >= 5e3 & x$KB < 10e3] * 1e3),
        sum(x$KB[x$KB >= 10e3] * 1e3)
    )

    return(d/autosome_len)
}



#ids <- read.table(fam, h=F, stringsAsFactors=F)$V1
ids <- c("9264_2", "dups1", "9259", "9253","9279", "9176merged", "P9", "P13merged")
idnew <- df$newPaperIDs
names(idnew) <- df$ID
ids2 <- idnew[ids]

pop <- c(rep("Northern", 5), rep("Southern",3))
saola_colors <- c("Northern" = "#86ABCB",
                  "Southern" = "#FFB300")
rohcolsnorth <- c("#bcd1e2", "#98b8d3", "#6392bc", "#385f82")
rohcolscent <- c("#ffe199", "#ffca4d", "#ffb300", "#b37d00")

rohcolslegend <- c("#cccccc", "#a6a6a6", "#8f8f8f", "#5d5d5d")

m <- sapply(ids, getROHproportion, roh=roh, autosome_len=autosome_len)
m2 <- rbind(m,m)
m2[1:4,6:8] <- 0
m2[5:8,1:5] <- 0

#
density <- c(0, 1, 2, 2)
angle <- c(0,0,0,45)

outpng <-  "/home/lpq293/mypopgen/saola/paperplots2/figure2/saola_rohs.png" # name of output png file

bitmap(outpng, h=5,w=5,res=300)
par(oma=c(3,1,0,11))
#par(oma=c(2,0,0,1.5))

x <- barplot(m2, names.arg=ids2, las=2,
     col=c(rohcolsnorth, rohcolscent), 
        ylab="Fraction of genome in ROH", main="",
         cex.lab=2.5, cex.main=1.8, xlim=c(0,9),
         space=c(0.1,0.1,0.1,0.1,0.1,1.,0.1,0.1),
         yaxt="n",cex.lab=2.3,cex.names=2.3,xpd=NA)#, density=density, angle=angle)

axis(2, at=seq(0, 0.4, by=0.1), xpd=NA,cex.axis=2)
segments(x0=-0.5, x1= 10, y0=seq(0, 0.4, by=0.1)[-1], lwd=0.3, xpd=T, lty=2)
#abline(h=seq(0, 0.4, by=0.1)[-1], lty=2, lwd=0.3)
barplot(m2, xlab="",  yaxt="n",
        ylab="", space=c(0.1,0.1,0.1,0.1,0.1,1.,0.1,0.1),
        names.arg=rep("", ncol(m2)),
        col=c(rohcolsnorth, rohcolscent), add=T)

ybot <- seq(0.1,0.2,length.out=4)
ytop <- seq(0.1,0.2,length.out=4) + 0.1/3
ymed <- (ybot + ytop) / 2
rect(xleft=10.5, xright=11.5, ybottom=ybot,
     ytop=ytop,
     col=rohcolslegend,xpd=NA)
text(x=13, y=ymed, labels=c("1-2.5 Mbp", "2.5-5Mbp", "5-10 Mbp", ">10 Mbp"), xpd=NA,cex=1.8)
text(font=2,x=12,y=ytop[4] + 0.02, labels="ROH length", xpd=NA,cex=2)


#barplot(m, names.arg=ids, las=2, col=RColorBrewer::brewer.pal("Set2", n=6), ylab="Genome fraction", main="Runs of homozygosity distribution", cex.lab=1.2, cex.main=1.4,space=c(0.1,0.1,0.1,0.1,0.1,1.8,0.1,0.1))

         
#text(y=-0.08, x=c(1.5,6), labels=c("Central", "North"), xpd=NA, cex=1.5)
dev.off()


outpdf <-  "/home/lpq293/mypopgen/saola/paperplots2/figure2/saola_rohs.pdf" # name of output png file

pdf(outpdf, h=10,w=10)
par(oma=c(3,1,0,13))
#par(oma=c(2,0,0,1.5))

x <- barplot(m2, names.arg=ids2, 
        las=2, col=c(rohcolsnorth, rohcolscent), 
        ylab="Fraction of genome in ROH", main="", 
        cex.lab=3, cex.main=1.4, xlim=c(0,9),
        space=c(0.1,0.1,0.1,0.1,0.1,1.,0.1,0.1),yaxt="n",
        cex.axis=2.3,cex.names=2.3,xpd=NA)#, density=density, angle=angle)

axis(2, at=seq(0, 0.4, by=0.1), xpd=NA,cex.axis=2.4)
segments(x0=-0.5, x1= 10, y0=seq(0, 0.4, by=0.1)[-1], lwd=0.3, xpd=T, lty=2)
#abline(h=seq(0, 0.4, by=0.1)[-1], lty=2, lwd=0.3)
barplot(m2, xlab="",  yaxt="n",
        ylab="", space=c(0.1,0.1,0.1,0.1,0.1,1.,0.1,0.1),
        names.arg=rep("", ncol(m2)),
        col=c(rohcolsnorth, rohcolscent), add=T)

ybot <- seq(0.1,0.2,length.out=4)
ytop <- seq(0.1,0.2,length.out=4) + 0.1/3
ymed <- (ybot + ytop) / 2
rect(xleft=10.5, xright=11.5, ybottom=ybot,
     ytop=ytop,
     col=rohcolslegend,xpd=NA)
text(x=13, y=ymed, labels=c("1-2.5 Mbp", "2.5-5Mbp", "5-10 Mbp", ">10 Mbp"), xpd=NA,cex=1.8)
text(font=2,x=12,y=ytop[4] + 0.02, labels="ROH length", xpd=NA,cex=2)


#barplot(m, names.arg=ids, las=2, col=RColorBrewer::brewer.pal("Set2", n=6), ylab="Genome fraction", main="Runs of homozygosity distribution", cex.lab=1.2, cex.main=1.4,space=c(0.1,0.1,0.1,0.1,0.1,1.8,0.1,0.1))

         
#text(y=-0.08, x=c(1.5,6), labels=c("Central", "North"), xpd=NA, cex=1.5)
dev.off()
