source("/home/lpq293/mypopgen/saola/info/loadPopInfo.R")

alphaCol <- function(col, transPct){

    rgbCol <- col2rgb(col)
    alphaCol <- rgb(rgbCol[1], rgbCol[2], rgbCol[3],
                    max = 255, alpha = (100 - transPct) * 255 / 100)
    return(alphaCol)
}



hetfile <- "/home/lpq293/mypopgen/saola/analyses/heterozygosities/collectplot/saola_hets_sfs.collected.txt"
hets <- read.table(hetfile, h=T, stringsAsFactors=F)
rownames(hets) <- hets$ID


### read in noroh plots
f1 <- "/home/lpq293/mypopgen//saola/analyses/heterozygosities/hets_noroh/results_allmuts/sfs/collected.txt"
f2 <- "/home/lpq293/mypopgen//saola/analyses/heterozygosities/hets_noroh/results_notrans/sfs/collected.txt"

hnoroh_allmuts <- read.table(f1, h=T, stringsAsFactors=F)
hnoroh_notrans <- read.table(f2, h=T, stringsAsFactors=F)

hets <- hets[df$ID[is.na(df$exclude_reason) & df$error_rates_transversions < 0.001],]

# merge hets with without roh to plot togetehr
h <- merge(hets,hnoroh_notrans[,-(2:4)], by.x="ID", by.y="id")
colnames(h) <- c(colnames(h)[-ncol(h)], "Het_notransitions_norohs")
# add pop info
names(pop) <- id
h$pop <- pop[h$ID]

ylim <- range(c(hets$Het_notransitions, h$Het_notransitions_norohs) * 3) * c(0.9999, 1.0001)


outpng <- "saola_hets.png"
#bitmap(outpng, w=5.5,h=5,res=300,type="pngalpha")
png(outpng,w=600,h=550,units="px")#,res=300)
par(oma=c(0,6,0,0))

boxplot(hets$Het_notransitions * 3 ~ factor(pop[hets$ID], levels=popord),
        ylim=ylim, ,yaxt="n", cex.axis=3,
        col=sapply(saola_colors, alphaCol, transPct=50), outline=F, xaxt="n",
        xlim=c(0,6), at=c(1,4), xlab="", ylab="", border=alphaCol("black", 50))

points(hets$Het_notransitions * 3 ~ jitter(ifelse(as.integer(factor(pop[hets$ID], levels=popord))==1,1,4),0.5),
       bg=sapply(saola_colors, alphaCol, transPct=70)[pop[hets$ID]], pch=c(21, 23)[hets$ID %in% h$ID + 1],
       cex=3, col=alphaCol("black", 50))

boxplot(h$Het_notransitions_norohs* 3 ~ factor(pop[h$ID], levels=popord),
        ylim=ylim, ,yaxt="n", cex.axis=3,
        col=sapply(saola_colors, alphaCol, transPct=0), outline=F, xaxt="n",
        xlim=c(0,6), at=c(2,5),add=T,xlab="",ylab="", border=alphaCol("black", 0))

points(h$Het_notransitions_norohs* 3 ~ jitter(ifelse(as.integer(factor(pop[h$ID], levels=popord))==1,2,5),0.5),
       bg=sapply(saola_colors, alphaCol, transPct=0)[pop[h$ID]], pch=23,
       cex=3, col=alphaCol("black", 0))


axis(2,labels=format(signif(seq(ylim[1],ylim[2], length.out=5),1),  scientific=FALSE),
     at = signif(seq(ylim[1],ylim[2], length.out=5),1), las=2, cex.axis=2)
title(ylab="Genome-wide heterozygozity",line=7, cex.lab=2.5,xpd=NA)

axis(1, at=c(1.5, 4.5), labels=FALSE, cex.axis=2.5)
axis(1, at=c(1.5, 4.5), labels=popord, tick=F, line=1.5, cex.axis=2.5)

# fake boxplots for legend
boxplot(c(0.000545, 0.00056, 0.000575) - 0.00001, at=2.92,add=T,xaxt="n", yaxt="n",
        col=alphaCol("darkgrey", 0), border=alphaCol("black", 0), cex=2)
boxplot(c(0.000545, 0.00056, 0.000575) + 0.00006, at=2.92, add=T,xaxt="n", yaxt="n",
        col=alphaCol("darkgrey", 70), border=alphaCol("black", 50), cex=2)
text(x=4.8,y=c( 0.00055,  0.00062),
     labels=c("Excluding ROHs",
              "Including ROHs"),
              cex=2)

dev.off()





outpdf <- "saola_hets.pdf"

#bitmap(outpng, w=5.5,h=5,res=300,type="pngalpha")
pdf(outpdf,w=8,h=7)#,res=300)
par(oma=c(0,6,0,0))

boxplot(hets$Het_notransitions * 3 ~ factor(pop[hets$ID], levels=popord),
        ylim=ylim, ,yaxt="n", cex.axis=3,
        col=sapply(saola_colors, alphaCol, transPct=50), outline=F, xaxt="n",
        xlim=c(0,6), at=c(1,4), xlab="", ylab="", border=alphaCol("black", 50))

points(hets$Het_notransitions * 3 ~ jitter(ifelse(as.integer(factor(pop[hets$ID], levels=popord))==1,1,4),0.5),
       bg=sapply(saola_colors, alphaCol, transPct=70)[pop[hets$ID]], pch=c(21, 23)[hets$ID %in% h$ID + 1],
       cex=3, col=alphaCol("black", 50))

boxplot(h$Het_notransitions_norohs* 3 ~ factor(pop[h$ID], levels=popord),
        ylim=ylim, ,yaxt="n", cex.axis=3,
        col=sapply(saola_colors, alphaCol, transPct=0), outline=F, xaxt="n",
        xlim=c(0,6), at=c(2,5),add=T,xlab="",ylab="", border=alphaCol("black", 0))

points(h$Het_notransitions_norohs* 3 ~ jitter(ifelse(as.integer(factor(pop[h$ID], levels=popord))==1,2,5),0.5),
       bg=sapply(saola_colors, alphaCol, transPct=0)[pop[h$ID]], pch=23,
       cex=3, col=alphaCol("black", 0))


axis(2,labels=format(signif(seq(ylim[1],ylim[2], length.out=5),1),  scientific=FALSE),
     at = signif(seq(ylim[1],ylim[2], length.out=5),1), las=2, cex.axis=2)
title(ylab="Genome-wide heterozygozity",line=7, cex.lab=2.5,xpd=NA)

axis(1, at=c(1.5, 4.5), labels=FALSE, cex.axis=2.5)
axis(1, at=c(1.5, 4.5), labels=popord, tick=F, line=1.5, cex.axis=2.5)

# fake boxplots for legend
boxplot(c(0.000545, 0.00056, 0.000575) - 0.00001, at=2.92,add=T,xaxt="n", yaxt="n",
        col=alphaCol("darkgrey", 0), border=alphaCol("black", 0), cex=2)
boxplot(c(0.000545, 0.00056, 0.000575) + 0.00006, at=2.92, add=T,xaxt="n", yaxt="n",
        col=alphaCol("darkgrey", 70), border=alphaCol("black", 50), cex=2)
text(x=4.8,y=c( 0.00055,  0.00062),
     labels=c("Excluding ROHs",
              "Including ROHs"),
              cex=2)

dev.off()

