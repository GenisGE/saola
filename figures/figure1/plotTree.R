library(ape)
library(strap)

f1 <- "/home/lpq293/mypopgen/saola/paperplots2/figure1/files/tree.nwk"
f2 <- "/home/lpq293/mypopgen/saola/paperplots2/figure1/files/kuduAge.txt"

#
sps <-  c("Minke Whale", "Goat", "Lesser Kudu", "Saola", "Water buffalo", "Cattle")

tree <- read.tree(f1)
times <- read.table(f2, h=T, as.is=T, row.names=1)
rownames(times) <- sps


sps <-  c("Minke Whale", "Goat", "Cattle", "Water buffalo", "Saola", "Lesser Kudu")

tree$tip.label <- sps

tree$root.time <- max(tree$edge.length)#max(times[,2:3])# * 1.05

times <- times[tree$tip.label,]

outpng <- "/home/lpq293/mypopgen/saola/paperplots2/figure1/tree.png"

bitmap(outpng, w=12, h=4, res=300)
par(mar=c(0,0,0,0))
geoscalePhylo(tree, 
    units=c("Period", "Epoch"), #units=c("Period", "Epoch"), 
    cex.age=3, cex.tip=2.5, boxes="Epoch", 
    cex.ts=3, ts.col=T, quat.rm=T, erotate=0, 
    width=4, font=2, label.offset=1, x.lim=c(-8,71))

offset <- max(tree$edge.length)

nodey <- c("Minke Whale" = 2.3125,
           "Goat" = 3.5625,
           "Cattle" = 3.5,
           "Saola" = 4.25,
           "Lesser Kudu" = 5.125)

do <- c("Minke Whale", "Goat", "Cattle", "Saola", "Lesser Kudu")
segments(x0= offset - times[do,]$LAD, x1=offset - times[do,]$FAD, y0=nodey[do], lwd=10)
dev.off()


outpdf <- "/home/lpq293/mypopgen/saola/paperplots2/figure1/tree.pdf"
pdf(outpdf, w=24, h=8)
#par(mar=c(0,0,0,0))
par(oma=c(1,0,0,0), xpd=NA)
geoscalePhylo(tree, 
    units=c("Period", "Epoch"), #units=c("Period", "Epoch"), 
    cex.age=3, cex.tip=2.5, boxes="Epoch", 
    cex.ts=3, ts.col=T, quat.rm=T, erotate=0, 
    width=4, font=2, label.offset=1, x.lim=c(-8,71))

offset <- max(tree$edge.length)

nodey <- c("Minke Whale" = 2.3125,
           "Goat" = 3.5625,
           "Cattle" = 3.5,
           "Saola" = 4.25,
           "Lesser Kudu" = 5.125)

do <- c("Minke Whale", "Goat", "Cattle", "Saola", "Lesser Kudu")
segments(x0= offset - times[do,]$LAD, x1=offset - times[do,]$FAD, y0=nodey[do], lwd=10)
dev.off()



if(FALSE){
### THIS WILL LOAD DATA FRAME WITH GEOLOGICAL TIMES CALE,
### FOR IN THE FUTURE MAKING IT MYSELF
vers = "ICS2013"
utils::data(timescales, envir = environment())
timescale <- timescales[[vers]]
}
