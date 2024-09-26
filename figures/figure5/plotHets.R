

saolaf <- "/home/lpq293/mypopgen/saola/analyses/heterozygosities/collectplot/saola_hets_sfs.collected.txt"
otherf <- "/home/lpq293/mypopgen/saola/external_resource/fromrhinocell/otherhets.tsv.csv"

saola_colors <- c( "Northern" = "#86ABCB",  "Southern" = "#FFB300")
#saola_colors <- c( "Northern" = "#29455e",  "Central" = "#996b00")


hets <- read.table(saolaf, h=T, stringsAsFactors=F)
rownames(hets) <- hets$ID

k <- c("9176merged", "P13merged", "P9", "9264_2", "dups1", "9259", "9253", "9279")

saolahets <- data.frame(id=k,
                        pop=c(rep("Southern", 3), rep("Northern", 5)),
                        het=hets[k, "Het_notransitions"] * 3,
                        stringsAsFactors=F)


otherhets <- read.delim(otherf, skip=1, h=1, stringsAsFactors=F, sep="\t")


mammal_orders <- c("Pilosa", "Carnivora", "Perissodactyla", "Artiodactyla",  "Rodentia", "Primates", "Chiroptera", "Eulipotyphla", "Lagomorpha", "Dermoptera", "Scandentia", "Cingulata", "Hyracoidea", "afrosoricida", "Macroscelidea", "Tubulidentata", "Dasyuromorphia")



outpng <- "hetsContext.png"
bitmap(outpng, w=10, h=3, res=300)

par(mfrow=c(1,2), oma = c(0,1,0,1))

##########################################
#### HISTOGRAM 1 COLORED BY TAXONOMY #####
##########################################
colpalother <- c("lightgrey", "firebrick3", "dodgerblue3")
colpalother <- c("lightgrey", "seagreen3", "sienna3")

cexlab <- 2.5
cexlegend <- 2
cexaxis <- 2
cexpt <- 2.5
hist(log10(otherhets$Heterozygosity),breaks=75,
     xlab="Genome-wide heterozygosity", xaxt="n", col=colpalother[1], border=NA, main="",
     cex.lab=cexlab, cex.axis=cexaxis, xpd=NA)

k <- otherhets$Family %in% mammal_orders

hist(log10(otherhets$Heterozygosity[k]),
     breaks=75, , xaxt="n", col=colpalother[2], border=NA, add=T)

k <- otherhets$Order == "Bovidae"

hist(log10(otherhets$Heterozygosity[k]),
     breaks=75, , xaxt="n", col=colpalother[3], border=NA, add=T)
hetaxis <- c(1e-4, 1e-3, 1e-2, 1e-1)
axis(1, at=log10(hetaxis), labels=format(hetaxis), cex.axis=cexaxis)

# mark position of some relevant species in plot
label1 <- c("Baiji", "Nilgiri Tahr", "Hirola")
aty = c(5,10,7)
atx = log10(c(1.2e-4, 2e-4, 4.2e-4))
text(y=aty, x = atx, labels=label1, cex=1.8)
arrows(x0=atx, 
     x1=log10(otherhets$Heterozygosity[sapply(label1, function(x) which(otherhets$Common.name == x))]),
     y0=aty - 0.75, y1=1,
     angle=20, length=0.07, lwd=1)

points(y=jitter(rep(-1, nrow(saolahets)),20),x=rev(log10(saolahets$het)), bg=rev(saola_colors[saolahets$pop]), cex=cexpt, pch=21, xpd=NA)

legend(y=25, x=log10(1.5e-2), legend=c("Animals", "Mammals", "Bovids"), title="Taxonomy", 
     fill=colpalother, bty="n", cex=cexlegend, xpd=NA)

legend(y=13, x=log10(1.5e-2), legend=names(saola_colors), pt.bg=saola_colors, bty="n",
       cex=cexlegend,
       title="Saola\nheterozygosity", pch=21, pt.cex=cexpt, xpd=NA)


#################################################
#### HISTOGRAM 2 COLORED BY IUCN ASSESSMENT #####
#################################################
colpalother <- c("lightgrey", "firebrick3", "dodgerblue3")
colpalother <- c("lightgrey", "indianred1", "indianred3", "indianred4")
colpalother <- c("lightgrey", "goldenrod2", "darkorange1", "firebrick4")


k <- otherhets$Conservation.status %in% c("Critically Endangered", "Endangered", "Vulnerable", "Near Threatened", "Least Concern")

hist(log10(otherhets$Heterozygosity[k]),breaks=100,
     xlab="Genome-wide heterozygosity", xaxt="n", col=colpalother[1], border=NA, main="",
     cex.lab=cexlab, cex.axis=cexaxis, xpd=NA)


k <- otherhets$Conservation.status %in% c("Critically Endangered", "Endangered", "Vulnerable")
#k <- otherhets$Conservation.status %in% c("Vulnerable")


hist(log10(otherhets$Heterozygosity[k]),
     breaks=100, , xaxt="n", col=colpalother[2], border=NA, add=T)
hetaxis <- c(1e-4, 1e-3, 1e-2, 1e-1)


k <- otherhets$Conservation.status %in% c("Critically Endangered", "Endangered")
#k <- otherhets$Conservation.status %in% c("Endangered")

hist(log10(otherhets$Heterozygosity[k]),
     breaks=100, , xaxt="n", col=colpalother[3], border=NA, add=T)


k <- otherhets$Conservation.status %in% c("Critically Endangered")

hist(log10(otherhets$Heterozygosity[k]),
     breaks=100, , xaxt="n", col=colpalother[4], border=NA, add=T)


axis(1, at=log10(hetaxis), labels=format(hetaxis), cex.axis=cexaxis)

# mark position of some relevant species in plot
label1 <- c("Cheetah", "Snow leopard", "Tasmanian devil")
aty = c(3.5,5.3,7)
atx = log10(c(1.42e-4, 2.2e-4, 5.1e-4))
text(y=aty, x = atx, labels=label1, cex=1.8)
arrows(x0=atx, 
     x1=log10(otherhets$Heterozygosity[sapply(label1, function(x) which(otherhets$Common.name == x))]),
     y0=aty - 0.75, y1=1,
     angle=20, length=0.07, lwd=1)


points(y=jitter(rep(-0.35, nrow(saolahets)),20),x=rev(log10(saolahets$het)), bg=rev(saola_colors[saolahets$pop]), cex=cexpt, pch=21, xpd=NA)

legend(y=11, x=log10(0.8e-2),
       legend=c("All with\nassessment",
                "Vulnerable",
                "Endangered",
                "Critically\nEndangered"),
       fill=colpalother, bty="n", cex=cexlegend, title="IUCN assessment", xpd=NA)



dev.off()








outpdf <- "hetsContext.pdf"
pdf(outpdf, w=20, h=6)

par(mfrow=c(1,2), oma = c(0,1,0,1))

##########################################
#### HISTOGRAM 1 COLORED BY TAXONOMY #####
##########################################
colpalother <- c("lightgrey", "firebrick3", "dodgerblue3")
colpalother <- c("lightgrey", "seagreen3", "sienna3")

cexlab <- 2.5
cexlegend <- 2
cexaxis <- 2
cexpt <- 2.5
hist(log10(otherhets$Heterozygosity),breaks=75,
     xlab="Genome-wide heterozygosity", xaxt="n", col=colpalother[1], border=NA, main="",
     cex.lab=cexlab, cex.axis=cexaxis, xpd=NA)

k <- otherhets$Family %in% mammal_orders

hist(log10(otherhets$Heterozygosity[k]),
     breaks=75, , xaxt="n", col=colpalother[2], border=NA, add=T)

k <- otherhets$Order == "Bovidae"

hist(log10(otherhets$Heterozygosity[k]),
     breaks=75, , xaxt="n", col=colpalother[3], border=NA, add=T)
hetaxis <- c(1e-4, 1e-3, 1e-2, 1e-1)
axis(1, at=log10(hetaxis), labels=format(hetaxis), cex.axis=cexaxis)

# mark position of some relevant species in plot
label1 <- c("Baiji", "Nilgiri Tahr", "Hirola")
aty = c(5,10,7)
atx = log10(c(1.2e-4, 2e-4, 4.2e-4))
text(y=aty, x = atx, labels=label1, cex=1.8)
arrows(x0=atx, 
     x1=log10(otherhets$Heterozygosity[sapply(label1, function(x) which(otherhets$Common.name == x))]),
     y0=aty - 0.75, y1=1,
     angle=20, length=0.07, lwd=1)

points(y=jitter(rep(-1, nrow(saolahets)),20),x=rev(log10(saolahets$het)), bg=rev(saola_colors[saolahets$pop]), cex=cexpt, pch=21, xpd=NA)

legend(y=25, x=log10(1.5e-2), legend=c("Animals", "Mammals", "Bovids"), title="Taxonomy", 
     fill=colpalother, bty="n", cex=cexlegend, xpd=NA)

legend(y=13, x=log10(1.5e-2), legend=names(saola_colors), pt.bg=saola_colors, bty="n",
       cex=cexlegend,
       title="Saola\nheterozygosity", pch=21, pt.cex=cexpt, xpd=NA)


#################################################
#### HISTOGRAM 2 COLORED BY IUCN ASSESSMENT #####
#################################################
colpalother <- c("lightgrey", "firebrick3", "dodgerblue3")
colpalother <- c("lightgrey", "indianred1", "indianred3", "indianred4")
colpalother <- c("lightgrey", "goldenrod2", "darkorange1", "firebrick4")


k <- otherhets$Conservation.status %in% c("Critically Endangered", "Endangered", "Vulnerable", "Near Threatened", "Least Concern")

hist(log10(otherhets$Heterozygosity[k]),breaks=100,
     xlab="Genome-wide heterozygosity", xaxt="n", col=colpalother[1], border=NA, main="",
     cex.lab=cexlab, cex.axis=cexaxis, xpd=NA)


k <- otherhets$Conservation.status %in% c("Critically Endangered", "Endangered", "Vulnerable")
#k <- otherhets$Conservation.status %in% c("Vulnerable")


hist(log10(otherhets$Heterozygosity[k]),
     breaks=100, , xaxt="n", col=colpalother[2], border=NA, add=T)
hetaxis <- c(1e-4, 1e-3, 1e-2, 1e-1)


k <- otherhets$Conservation.status %in% c("Critically Endangered", "Endangered")
#k <- otherhets$Conservation.status %in% c("Endangered")

hist(log10(otherhets$Heterozygosity[k]),
     breaks=100, , xaxt="n", col=colpalother[3], border=NA, add=T)


k <- otherhets$Conservation.status %in% c("Critically Endangered")

hist(log10(otherhets$Heterozygosity[k]),
     breaks=100, , xaxt="n", col=colpalother[4], border=NA, add=T)


axis(1, at=log10(hetaxis), labels=format(hetaxis), cex.axis=cexaxis)

# mark position of some relevant species in plot
label1 <- c("Cheetah", "Snow leopard", "Tasmanian devil")
aty = c(3.5,5.3,7)
atx = log10(c(1.42e-4, 2.2e-4, 5.1e-4))
text(y=aty, x = atx, labels=label1, cex=1.8)
arrows(x0=atx, 
     x1=log10(otherhets$Heterozygosity[sapply(label1, function(x) which(otherhets$Common.name == x))]),
     y0=aty - 0.75, y1=1,
     angle=20, length=0.07, lwd=1)


points(y=jitter(rep(-0.35, nrow(saolahets)),20),x=rev(log10(saolahets$het)), bg=rev(saola_colors[saolahets$pop]), cex=cexpt, pch=21, xpd=NA)

legend(y=11, x=log10(0.8e-2),
       legend=c("All with\nassessment",
                "Vulnerable",
                "Endangered",
                "Critically\nEndangered"),
       fill=colpalother, bty="n", cex=cexlegend, title="IUCN assessment", xpd=NA)



dev.off()
