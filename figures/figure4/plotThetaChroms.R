source("/home/lpq293/mypopgen/saola/info/loadPopInfo.R")

read_merge_tables <- function(inprefix, minqsites_keep = 0.05, gfile =  "/home/lpq293/mypopgen/saola/genome_masks/v2/saola_ref/all/beds/saolaRefN_autosomal_10MBscaff.bed"){

    f1 <- paste0(inprefix, "_pi.txt")
    f2 <- paste0(inprefix,"_dxy.txt")
    f3 <- paste0(inprefix, "_combinedPops_pi.txt")


    ## read in data
    pi <- read.table(f1, h=T, stringsAsFactors=F)

    dxy <- read.table(f2, h=T, stringsAsFactors=F)

    picomb <- read.table(f3, h=T, stringsAsFactors=F)


    # if there is genome file, keep only chromosomes present there
    if(!is.null(gfile)){
        g <- read.table(gfile, h=F, stringsAsFactors=F)
        chroms <- g$V1
        pi <- pi[pi$chromosome %in% chroms,]
        dxy <- dxy[dxy$chromosome %in% chroms,]
        picomb <- picomb[picomb$chromosome %in% chroms,]
    }

    # find quantile of number of sites to remove
    minsitespi <- quantile(pi$no_sites, minqsites_keep)
    minsitesdxy <- quantile(dxy$no_sites, minqsites_keep)
    minsitespi2 <- quantile(picomb$no_sites, minqsites_keep)
    # clean out windows with too low number of sites
    pi <- pi[pi$no_sites > minsitespi,]         
    dxy <- dxy[dxy$no_sites > minsitesdxy,]         
    picomb <- picomb[picomb$no_sites > minsitespi2,]

    ### put together all data in single data frame
    pi1 <- pi[pi$pop == "Northern", ]
    pi2 <- pi[pi$pop == "Central", ]

    pinew <- merge(pi1, pi2, by=c("chromosome", "window_pos_1", "window_pos_2"), suffixes=c(".northern", ".central"))
    d <- merge(pinew, dxy, by=c("chromosome", "window_pos_1", "window_pos_2"), suffixes=c("", ".dxy"))
    d <- merge(d, picomb, by=c("chromosome", "window_pos_1", "window_pos_2"), suffixes=c(".dxy", ".combined"))
    ### resort data because merge functions messed order
    d <- d[order(d$chromosome, d$window_pos_1),]
    d <- d[unlist(sapply(chroms, function(x) which(d$chromosome == x))),]

    return(d)

    
}


saola_colors <- c(saola_colors, "Combined" = "#505050")
saola_colors2 <- c("Northern" = "#325675", "Central" = "#7f5900", "Combined" = "black")
saola_colors2 <- c("Northern" = "#325675", "Southern" = "#7f5900", "Combined" = "black")

northpal <- colorRampPalette(c("white", "#86ABCB"))
centpal <- colorRampPalette(c("white", "#FFB300"))
allpal <- colorRampPalette(c("white", "#505050"))

gfile <-  "/home/lpq293/mypopgen/saola/genome_masks/v2/saola_ref/all/beds/saolaRefN_autosomal_10MBscaff.bed"
inprefix <- "/home/lpq293/mypopgen/saola/analyses/theta_win/results_matchedsamples/pixy_out/output_maxmiss0.5_win100000"
outpre <- "/home/lpq293/mypopgen/saola/paperplots2/figure4/allscaffthetaplots/thetaplot_miss05_win100000"

rohfile1 <- "/home/lpq293/mypopgen/saola/paperplots/figure4/rohintersections/roh_intersect_northern.bed"
rohfile2 <- "/home/lpq293/mypopgen/saola/paperplots/figure4/rohintersections/roh_intersect_central.bed"
rohfile3 <- "/home/lpq293/mypopgen/saola/paperplots/figure4/rohintersections/roh_intersect_all.bed"

### READ IN DATA ####
rohsnorth <- read.table(rohfile1, h=F, stringsAsFactors=F)
rohscentral <- read.table(rohfile2, h=F, stringsAsFactors=F)
rohsall <- read.table(rohfile3, h=F, stringsAsFactors=F)

# sample size for each pop
n_north <- 3
n_central <- 3
n_all <- n_north + n_central

g <- read.table(gfile, h=F, stringsAsFactors=F)

d <- read_merge_tables(inprefix)
d$midpos <- (d$window_pos_1 +  d $window_pos_2)  / 2

#quantile(c(d$avg_pi.northern, d$avg_pi.central, d$avg_pi), 0.999)
ylims <- range(c(d$avg_pi.northern, d$avg_pi.central, d$avg_pi))

# function to generate histogram data without for vector x, restircted to ylims
# values above ylims are set to max lim
doHistLims <- function(x, ylims, nbreaks=100, rmout=TRUE){

    hsbreaks <- seq(0, ylims[2], length.out=nbreaks)
    if(!rmout) x[x > ylims[2]] <- ylims[2]
    if(rmout) x <- x[x < ylims[2]]
    hs <- hist(x, breaks=hsbreaks, plot=FALSE)
    return(hs)

}



hsbreaks <- seq(0, max(d$avg_pi), length.out=100)

chroms <- unique(d$chromosome)

c <- "1b"


#for(c in chroms){
outpng <- paste0("divacrosscrhrom_2b.png")
###outpng <- "test.png"
scaffsize <- c(g[g$V1==c,2], g[g$V1==c,3])
k <- d$chromosome == c

# get common ylims anc position for yaxt
ylims <- c(0, max(c(d$avg_pi.northern[k], d$avg_pi.central[k], d$avg_pi[k])))
yaxtpos <- signif(seq(0, ylims[2], length.out=5),2)

# prepare histograms with genome-wide distribution
# force breaks in hist to be the same and consistent with the range of pi in the plotted scaffold
hsbreaks <- seq(0, ylims[2], length.out=100)

hsnorth <- doHistLims(d$avg_pi.northern, ylims)
hscentral <- doHistLims(d$avg_pi.central, ylims)
hscombined <- doHistLims(d$avg_pi, ylims)

maxhsfreq <- max(c(hsnorth$counts/sum(hsnorth$counts),
                    hscentral$counts/sum(hscentral$counts),
                    hscombined$counts/sum(hscombined$counts)))

cex.axis <- 2.75
cex.lab <- 3.5
cex.lab2 <- 3

bitmap(outpng, w=12,h=3,res=300)
par(mfrow=c(3,1), mar=c(0,5.5,2,8), oma=c(8,5,7,30))#, oma=c(0,0,8,0))

####### Plot northern populaiton pi########
plot(0, type="n",
     xlab="", ylab="", xaxt="n", yaxt="n", bty="n",
     ylim=ylims, xlim=scaffsize)#c(0,sum(k)))

## BUT FIRST PLOT ROH LEGEND ABOVE THE PLOT
xpos <- scaffsize[2] * c(0.8, 1)
ylegend <- ylims[2] * c(1.4, 1.8)
ypos <- seq(ylegend[1], ylegend[2], length.out=4)

rect(ybottom=ypos[3], ytop=ypos[4],
     xleft = seq(xpos[1], xpos[2], length.out=5)[-5],
     xright = seq(xpos[1], xpos[2], length.out=5)[-1],
     col=northpal(n_north + 1), xpd=NA, border=NA)
rect(ybottom=ypos[2], ytop=ypos[3],
     xleft = seq(xpos[1], xpos[2], length.out=5)[-5],
     xright = seq(xpos[1], xpos[2], length.out=5)[-1],
     col=centpal(n_central + 1), xpd=NA, border=NA)
rect(ybottom=ypos[1], ytop=ypos[2],
     xleft = seq(xpos[1], xpos[2], length.out=8)[-8],
     xright = seq(xpos[1], xpos[2], length.out=8)[-1],
     col=allpal(n_all + 1), xpd=NA, border=NA)
rect(ybottom=ylegend[1], ytop=ylegend[2],
     xleft = xpos[1],
     xright = xpos[2],
     xpd=NA, lwd=0.5)
text(x=c(xpos[1], mean(xpos), xpos[2]), y=ylims[2] * 1.25, labels=c(0,0.5,1), cex=2.5, xpd=NA)
text(x=mean(xpos), y=ylims[2] * 1.95, labels="Local ROH density", cex=2.75, xpd=NA)

### NOT PLOT FOR REAL NORHTERN POPULATION PI####

#draw rohs
kroh <- rohsnorth$V1 == c
if(sum(kroh) > 0)
     rect(xleft=rohsnorth$V2[kroh], xright=rohsnorth$V3[kroh], 
          ybottom=0,ytop=ylims[2], border=0, 
          col=northpal(n_north+1)[rohsnorth$V4[kroh] + 1])

# draw points
points(y=d$avg_pi.northern[k], x=d$midpos[k], col=saola_colors2["Northern"], type="b", lwd=1, pch=16)
# draw y axis
axis(2, at =yaxtpos,  labels=yaxtpos, las=2, line=-2, cex.axis=cex.axis)
title(ylab=expression(Pi[Northern]),
           cex.lab=cex.lab2, line=7.5, xpd=NA)

## draw histogram
hs <- hsnorth
segments(x0=scaffsize[2] * 1.03,
          x1=scaffsize[2] * 1.03 + hs$counts/sum(hs$counts) * scaffsize[2], # rescale hist count values
          y0=hs$mids, y1=hs$mids,
          lwd=2, col=saola_colors["Northern"], xpd=NA)


####### Plot central populaiton pi########
plot(0, type="n",
     xlab="", ylab="", xaxt="n", yaxt="n", bty="n",
     ylim=ylims, xlim=scaffsize)#c(0,sum(k)))
#draw rohs
kroh <- rohscentral$V1 == c
if(sum(kroh) > 0)
     rect(xleft=rohscentral$V2[kroh], xright=rohscentral$V3[kroh], ybottom=0,ytop=ylims[2], border=0, col=centpal(n_central+1)[rohscentral$V4[kroh] + 1])

# draw points
points(y=d$avg_pi.central[k], x=d$midpos[k], col=saola_colors2["Southern"], type="b", lwd=1, pch=16)
# draw y axis
axis(2, at =yaxtpos,  labels=yaxtpos, las=2, line=-2, cex.axis=cex.axis)
title(ylab=expression(Pi[Southern]),
           cex.lab=cex.lab2, line=7.5, xpd=NA)

## draw histogram
hs <- hscentral
segments(x0=scaffsize[2] * 1.03,
          x1=scaffsize[2] * 1.03 + hs$counts/sum(hs$counts) * scaffsize[2], # rescale hist count values
          y0=hs$mids, y1=hs$mids,
          lwd=2, col=saola_colors["Southern"], xpd=NA)


#### Plot combined populations pi ########
plot(0, type="n",
     xlab="", ylab="", xaxt="n", yaxt="n", bty="n",
     ylim=ylims, xlim=scaffsize)#c(0,sum(k)))
#draw rohs
kroh <- rohsall$V1 == c
if(sum(kroh) > 0)
     rect(xleft=rohsall$V2[kroh], xright=rohsall$V3[kroh], ybottom=0,ytop=ylims[2], border=0, col=allpal(n_all+1)[rohsall$V4[kroh] + 1])


# draw points
points(y=d$avg_pi[k], x=d$midpos[k], col=saola_colors2["Combined"], type="b", lwd=1, pch=16)
# draw y axis
axis(2, at =yaxtpos,  labels=yaxtpos, las=2, line=-2, cex.axis=cex.axis)
title(ylab=expression(Pi[Combined]), cex.lab=cex.lab2, 
          line=7.5, xpd=NA)


## draw histogram
hs <- hscombined
segments(x0=scaffsize[2] * 1.03,
          x1=scaffsize[2] * 1.03 + hs$counts/sum(hs$counts) * scaffsize[2], # rescale hist count values
          y0=hs$mids, y1=hs$mids,
          lwd=2, col=saola_colors["Combined"], xpd=NA)



### draw x axes
## draw histogram x axis
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)

atlabs <- seq(0, ceiling_dec(maxhsfreq, 2), length.out=5)
atpos <- scaffsize[2] * 1.03 + atlabs * scaffsize[2]

axis(1, at = atpos, labels=FALSE, xpd=NA, line=1, cex.axis=cex.axis)
axis(1, at = atpos, tick=FALSE, labels=round(atlabs,2), xpd=NA, line=2, cex.axis=cex.axis)
text(x=mean(atpos), y=min(hs$mids) - ylims[2] * 0.95, 
     labels="Genome-wide window frequency", xpd=NA, cex=cex.lab)


# x axis wfor main plot wiht scaffold positions
axis(side=1, at = seq(scaffsize[1], scaffsize[2], length.out=10),
     labels = FALSE,
     line=1, xpd=NA, cex.axis=cex.axis)
axis(side=1, at = seq(scaffsize[1], scaffsize[2], length.out=10),
     tick=FALSE,
     labels = round(seq(scaffsize[1], scaffsize[2], length.out=10)/1e6),
     line=2, xpd=NA, cex.axis=cex.axis)
title(xlab=paste("Position (Mbp) in Predicted Chromosome Fragment",c), 
               line=6, xpd=NA, cex.lab=cex.lab)



dev.off()








#for(c in chroms){
outpdf <- "divacrosscrhrom_2b.pdf"
###outpng <- "test.png"
scaffsize <- c(g[g$V1==c,2], g[g$V1==c,3])
k <- d$chromosome == c

# get common ylims anc position for yaxt
ylims <- c(0, max(c(d$avg_pi.northern[k], d$avg_pi.central[k], d$avg_pi[k])))
yaxtpos <- signif(seq(0, ylims[2], length.out=5),2)

# prepare histograms with genome-wide distribution
# force breaks in hist to be the same and consistent with the range of pi in the plotted scaffold
hsbreaks <- seq(0, ylims[2], length.out=100)

hsnorth <- doHistLims(d$avg_pi.northern, ylims)
hscentral <- doHistLims(d$avg_pi.central, ylims)
hscombined <- doHistLims(d$avg_pi, ylims)

maxhsfreq <- max(c(hsnorth$counts/sum(hsnorth$counts),
                    hscentral$counts/sum(hscentral$counts),
                    hscombined$counts/sum(hscombined$counts)))

cex.axis <- 2.75
cex.lab <- 3.5
cex.lab2 <- 3



pdf(outpdf, w=24,h=6)
par(mfrow=c(3,1), mar=c(0,5.5,2,8), oma=c(8,5,7,30))#, oma=c(0,0,8,0))

####### Plot northern populaiton pi########
plot(0, type="n",
     xlab="", ylab="", xaxt="n", yaxt="n", bty="n",
     ylim=ylims, xlim=scaffsize)#c(0,sum(k)))

## BUT FIRST PLOT ROH LEGEND ABOVE THE PLOT
xpos <- scaffsize[2] * c(0.8, 1)
ylegend <- ylims[2] * c(1.4, 1.8)
ypos <- seq(ylegend[1], ylegend[2], length.out=4)

rect(ybottom=ypos[3], ytop=ypos[4],
     xleft = seq(xpos[1], xpos[2], length.out=5)[-5],
     xright = seq(xpos[1], xpos[2], length.out=5)[-1],
     col=northpal(n_north + 1), xpd=NA, border=NA)
rect(ybottom=ypos[2], ytop=ypos[3],
     xleft = seq(xpos[1], xpos[2], length.out=5)[-5],
     xright = seq(xpos[1], xpos[2], length.out=5)[-1],
     col=centpal(n_central + 1), xpd=NA, border=NA)
rect(ybottom=ypos[1], ytop=ypos[2],
     xleft = seq(xpos[1], xpos[2], length.out=8)[-8],
     xright = seq(xpos[1], xpos[2], length.out=8)[-1],
     col=allpal(n_all + 1), xpd=NA, border=NA)
rect(ybottom=ylegend[1], ytop=ylegend[2],
     xleft = xpos[1],
     xright = xpos[2],
     xpd=NA, lwd=0.5)
text(x=c(xpos[1], mean(xpos), xpos[2]), y=ylims[2] * 1.25, labels=c(0,0.5,1), cex=2.5, xpd=NA)
text(x=mean(xpos), y=ylims[2] * 1.95, labels="Local ROH density", cex=2.75, xpd=NA)

### NOT PLOT FOR REAL NORHTERN POPULATION PI####

#draw rohs
kroh <- rohsnorth$V1 == c
if(sum(kroh) > 0)
     rect(xleft=rohsnorth$V2[kroh], xright=rohsnorth$V3[kroh], 
          ybottom=0,ytop=ylims[2], border=0, 
          col=northpal(n_north+1)[rohsnorth$V4[kroh] + 1])

# draw points
points(y=d$avg_pi.northern[k], x=d$midpos[k], col=saola_colors2["Northern"], type="b", lwd=1, pch=16)
# draw y axis
axis(2, at =yaxtpos,  labels=yaxtpos, las=2, line=-2, cex.axis=cex.axis)
title(ylab=expression(Pi[Northern]),
           cex.lab=cex.lab2, line=7.5, xpd=NA)

## draw histogram
hs <- hsnorth
segments(x0=scaffsize[2] * 1.03,
          x1=scaffsize[2] * 1.03 + hs$counts/sum(hs$counts) * scaffsize[2], # rescale hist count values
          y0=hs$mids, y1=hs$mids,
          lwd=2, col=saola_colors["Northern"], xpd=NA)


####### Plot central populaiton pi########
plot(0, type="n",
     xlab="", ylab="", xaxt="n", yaxt="n", bty="n",
     ylim=ylims, xlim=scaffsize)#c(0,sum(k)))
#draw rohs
kroh <- rohscentral$V1 == c
if(sum(kroh) > 0)
     rect(xleft=rohscentral$V2[kroh], xright=rohscentral$V3[kroh], ybottom=0,ytop=ylims[2], border=0, col=centpal(n_central+1)[rohscentral$V4[kroh] + 1])

# draw points
points(y=d$avg_pi.central[k], x=d$midpos[k], col=saola_colors2["Southern"], type="b", lwd=1, pch=16)
# draw y axis
axis(2, at =yaxtpos,  labels=yaxtpos, las=2, line=-2, cex.axis=cex.axis)
title(ylab=expression(Pi[Southern]),
           cex.lab=cex.lab2, line=7.5, xpd=NA)

## draw histogram
hs <- hscentral
segments(x0=scaffsize[2] * 1.03,
          x1=scaffsize[2] * 1.03 + hs$counts/sum(hs$counts) * scaffsize[2], # rescale hist count values
          y0=hs$mids, y1=hs$mids,
          lwd=2, col=saola_colors["Southern"], xpd=NA)


#### Plot combined populations pi ########
plot(0, type="n",
     xlab="", ylab="", xaxt="n", yaxt="n", bty="n",
     ylim=ylims, xlim=scaffsize)#c(0,sum(k)))
#draw rohs
kroh <- rohsall$V1 == c
if(sum(kroh) > 0)
     rect(xleft=rohsall$V2[kroh], xright=rohsall$V3[kroh], ybottom=0,ytop=ylims[2], border=0, col=allpal(n_all+1)[rohsall$V4[kroh] + 1])


# draw points
points(y=d$avg_pi[k], x=d$midpos[k], col=saola_colors2["Combined"], type="b", lwd=1, pch=16)
# draw y axis
axis(2, at =yaxtpos,  labels=yaxtpos, las=2, line=-2, cex.axis=cex.axis)
title(ylab=expression(Pi[Combined]), cex.lab=cex.lab2, 
          line=7.5, xpd=NA)


## draw histogram
hs <- hscombined
segments(x0=scaffsize[2] * 1.03,
          x1=scaffsize[2] * 1.03 + hs$counts/sum(hs$counts) * scaffsize[2], # rescale hist count values
          y0=hs$mids, y1=hs$mids,
          lwd=2, col=saola_colors["Combined"], xpd=NA)



### draw x axes
## draw histogram x axis
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)

atlabs <- seq(0, ceiling_dec(maxhsfreq, 2), length.out=5)
atpos <- scaffsize[2] * 1.03 + atlabs * scaffsize[2]

axis(1, at = atpos, labels=FALSE, xpd=NA, line=1, cex.axis=cex.axis)
axis(1, at = atpos, tick=FALSE, labels=round(atlabs,2), xpd=NA, line=2, cex.axis=cex.axis)
text(x=mean(atpos), y=min(hs$mids) - ylims[2] * 0.95, 
     labels="Genome-wide window frequency", xpd=NA, cex=cex.lab)


# x axis wfor main plot wiht scaffold positions
axis(side=1, at = seq(scaffsize[1], scaffsize[2], length.out=10),
     labels = FALSE,
     line=1, xpd=NA, cex.axis=cex.axis)
axis(side=1, at = seq(scaffsize[1], scaffsize[2], length.out=10),
     tick=FALSE,
     labels = round(seq(scaffsize[1], scaffsize[2], length.out=10)/1e6),
     line=2, xpd=NA, cex.axis=cex.axis)
title(xlab=paste("Position (Mbp) in Predicted Chromosome Fragment",c), 
               line=6, xpd=NA, cex.lab=cex.lab)



dev.off()

