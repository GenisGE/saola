

otherf1 <- "/home/lpq293/mypopgen/saola/external_resource/fromrhinocell/otherhets.tsv.csv"
otherf2 <- "/home/lpq293/mypopgen/saola/external_resource/fromrhinocell/othergenload.tsv.csv"

otherhet <- read.delim(otherf1, skip=1, h=1, stringsAsFactors=F, sep="\t")
otherload <- read.delim(otherf2, skip=1, h=1, stringsAsFactors=F, sep="\t")

otherload$missenseratio <- (otherload$X..Missense + otherload$X..LoF.mutation)/ otherload$X..Silent
otherload$conservation.status[otherload$Common.name == "Asian elephant"] <- "Endangered"

f <- "/home/lpq293/mypopgen/saola/analyses/genetic_load/results/annotation_load/variant_impact_count_set3a.tsv"
saolaload <- read.table(f, h=T, stringsAsFactors=F)
f <- "/home/lpq293/mypopgen/saola/analyses/heterozygosities/collectplot/saola_hets_sfs.collected.txt"
saolahet <- read.table(f, h=T, stringsAsFactors=F)

k <- c("9176merged", "9264_2", "dups1")
rownames(saolahet) <- saolahet$ID
saolahet <- saolahet[k,]
#saolahet <- saolahet[saolahet$ID %in% k,]

saolamiss <- (saolaload$total_het_sites[saolaload$category == "MODERATE"] + saolaload$total_het_sites[saolaload$category == "HIGH"])/  saolaload$total_het_sites[saolaload$category == "LOW"]

# load collected het + missense ratio. quickly done many stuff missing
a <- read.table("collectedquick.tsv", h=T, stringsAsFactors=F, sep="\t")
a <- a[!is.na(a$het),]


a$cons <- sapply(a$commonname, function(x) otherload$conservation.status[otherload$Common.name == x][1])
meanhetsaola <- mean(saolahet$Het_allmuts)
meanmissratsaola <- mean(saolamiss)
a$commonname <- stringr::str_to_title(a$commonname)


otherload$Common.name <- stringr::str_to_title(otherload$Common.name)
otherload$Common.name <- gsub("Nigeran-Cameroon", "Nig.Cam.", otherload$Common.name)
otherload$Common.name <- gsub("Northern White Cheeked", "N. Wh. Chk.", otherload$Common.name)
otherload$Common.name <- gsub("Himalayan Field", "Him. Field", otherload$Common.name)
otherload$Common.name <- gsub("Great Roundleaf", "Gr. Rdlf.", otherload$Common.name)
otherload$Common.name <- gsub("Central", "Cnt.", otherload$Common.name)
otherload$Common.name <- gsub("Western", "West.", otherload$Common.name)
otherload$Common.name <- gsub("Eastern", "East.", otherload$Common.name)
otherload$Common.name <- gsub("Sumatran", "Sutan.", otherload$Common.name)
otherload$Common.name <- gsub("Borneon", "Born.", otherload$Common.name)
otherload$Common.name <- gsub("Siamang", "Siang.", otherload$Common.name)
otherload$Common.name <- gsub("Chimpanzee", "Chimp.", otherload$Common.name)
otherload$Common.name <- gsub("Orangutan", "Orang.", otherload$Common.name)
otherload$Common.name <- gsub("Pileated", "Pil.", otherload$Common.name)
otherload$Common.name <- gsub("African", "Afr.", otherload$Common.name)
otherload$Common.name <- gsub("Gibbon", "Gibb.", otherload$Common.name)


sps <- unique(otherload$Common.name)
conss <- sapply(sps, function(x) otherload$conservation.status[otherload$Common.name==x][1])

otherload$conservation.status  <- gsub("Least concern", "Least Concern", otherload$conservation.status)
conss <- gsub("Least concern", "Least Concern", conss)
a$cons <- gsub("Least concern", "Least Concern", a$cons)

conss <- factor(conss, levels= c("Least Concern", "Near Threatened", "Vulnerable",
                                 "Endangered", "Critically Endangered"))


meanmissratio <- tapply(otherload$missenseratio, factor(otherload$Common.name, levels=sps), mean)
ord <- order(conss, meanmissratio)



saola_colors <- c("Northern" = "#86ABCB", "Central" =  "#FFB300" )
colpalother <- c("Least Concern" = "lightgrey",
                 "Near Threatened" = "lightgrey",
                 "Vulnerable" = "goldenrod2",
                 "Endangered" = "darkorange1",
                 "Critically Endangered" = "firebrick4")

outpng <- "loadratios.png"

cexlab <- 2.5
cexlegend <- 2
cexaxis <- 2
cexpt <- 3

bitmap(outpng, h=3,w=10,res=300)
ylim <- range(c(otherload$missenseratio,saolamiss))
layout(t(c(1,1,2)))
par(oma = c(9.5,7,0,0))
boxplot(c(otherload$missenseratio,saolamiss)  ~ factor(c(otherload$Common.name, rep("Saola", 3)), levels = c(sps[ord], "Saola")),
        outline=F, las=2, ylab="",
        cex.lab=cexlab, xpd=NA, cex.axis=cexaxis, xlab="")
title(ylab="# LoF and non-synonymous hets. to\n# synonymous hets.", xpd=NA, cex.lab=cexlab, line=4.5)
#        col=c(colpalother[conss[ord]], "black"),
#        las=2)
points(c(otherload$missenseratio,saolamiss)  ~ jitter(as.integer(factor(c(otherload$Common.name, rep("Saola", 3)), levels = c(sps[ord], "Saola"))),1), pch=21,
       bg=c(colpalother[otherload$conservation.status], saola_colors[c("Central", "Northern", "Northern")]), cex=cexpt)


legend(x=0.5,y=1.5,title="IUCN assessment",
       pch=21,
       pt.bg= c("lightgrey", "goldenrod2", "darkorange1", "firebrick4"),
       legend=c("Least Concern or Near Threatened", "Vulnerable", "Endangered", "Critically Endangered"), cex=2, pt.cex=cexpt, bty="n")

par(mar=c(4,1,4,2))


plot(1, type="n", ylim=ylim, xlim=range(a$het) * c(0.9, 1.1),
     #ylab="# Missense heterozygotes / \n# Silent heterozygotes",
     ylab="",
     xlab="",
     xpd=NA, cex.lab=cexlab, cex.axis=cexaxis)

#text(x=a$het, y=a$missratio, labels=a$commonname)
#points(c(a$missratio, meanmissratsaola) ~ c(a$het, meanhetsaola),
#       pch=21,
#       bg=c(colpalother[a$cons], "purple"), cex=cexpt)
points(a$missratio ~ a$het,
       pch=21,
       bg=c(colpalother[a$cons], "purple"), cex=cexpt)
points(meanmissratsaola ~ meanhetsaola, pch=8, bg="black", cex=cexpt, lwd=2)
text(x=meanhetsaola * 3.4, y=meanmissratsaola * 1.06, labels="Saola", cex=3)
title(xlab="Genome-wide heterozygosity", xpd=NA, cex.lab=3, cexlab, line=4)

dev.off()







outpdf <- "loadratios.pdf"

cexlab <- 2.5
cexlegend <- 2
cexaxis <- 2
cexpt <- 3

pdf(outpdf, h=6,w=20)
ylim <- range(c(otherload$missenseratio,saolamiss))
layout(t(c(1,1,2)))
par(oma = c(9.5,7,0,0))
boxplot(c(otherload$missenseratio,saolamiss)  ~ factor(c(otherload$Common.name, rep("Saola", 3)), levels = c(sps[ord], "Saola")),
        outline=F, las=2, ylab="",
        cex.lab=cexlab, xpd=NA, cex.axis=cexaxis, xlab="")
title(ylab="# LoF and non-synonymous hets. to\n# synonymous hets.", xpd=NA, cex.lab=cexlab, line=4.5)
#        col=c(colpalother[conss[ord]], "black"),
#        las=2)
points(c(otherload$missenseratio,saolamiss)  ~ jitter(as.integer(factor(c(otherload$Common.name, rep("Saola", 3)), levels = c(sps[ord], "Saola"))),1), pch=21,
       bg=c(colpalother[otherload$conservation.status], saola_colors[c("Central", "Northern", "Northern")]), cex=cexpt)


legend(x=0.5,y=1.5,title="IUCN assessment",
       pch=21,
       pt.bg= c("lightgrey", "goldenrod2", "darkorange1", "firebrick4"),
       legend=c("Least Concern or Near Threatened", "Vulnerable", "Endangered", "Critically Endangered"), cex=2, pt.cex=cexpt, bty="n")

par(mar=c(4,1,4,2))


plot(1, type="n", ylim=ylim, xlim=range(a$het) * c(0.9, 1.1),
     #ylab="# Missense heterozygotes / \n# Silent heterozygotes",
     ylab="",
     xlab="",
     xpd=NA, cex.lab=cexlab, cex.axis=cexaxis)

#text(x=a$het, y=a$missratio, labels=a$commonname)
#points(c(a$missratio, meanmissratsaola) ~ c(a$het, meanhetsaola),
#       pch=21,
#       bg=c(colpalother[a$cons], "purple"), cex=cexpt)
points(a$missratio ~ a$het,
       pch=21,
       bg=c(colpalother[a$cons], "purple"), cex=cexpt)
points(meanmissratsaola ~ meanhetsaola, pch=8, bg="black", cex=cexpt, lwd=2)
text(x=meanhetsaola * 3.4, y=meanmissratsaola * 1.06, labels="Saola", cex=3)
title(xlab="Genome-wide heterozygosity", xpd=NA, cex.lab=3, cexlab, line=4)

dev.off()
