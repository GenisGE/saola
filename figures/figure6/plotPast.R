source("/home/lpq293/mypopgen/saola/info/loadPopInfo.R")

library(colorspace)



nruns <- 200
infolder <- "/home/lpq293/mypopgen/saola/analyses/simulations/nwf_sims/longer/perezmodel/" # nolint: line_length_linter.

read1RunRes <- function(i, infolder=infolder)
    read.table(paste0(infolder, "/", i, "/load_decomposed_byh.txt"), h=F, skip=1, stringsAsFactors=F, sep=",")[,-(1:3)] # nolint: line_length_linter.


read1RunLoadH <- function(i, infolder=infolder)
    read.table(paste0(infolder, "/", i, "/load_decomposed_byh.txt"), h=F, skip=1, stringsAsFactors=F, sep=",")[,-(1:3)]


allres <- lapply(1:nruns, read1RunLoadH, infolder=infolder)

binStarts <- c(0, 0.01, 0.05, seq(0.1, 0.9, by=0.1))
binEnds <- c(0.01, 0.05, seq(0.1, 1, by=0.1))
nbins <- 12



cols1 <- read.table(paste0(infolder, "/", 1, "/load_decomposed_byh.txt"),
    h=F, skip=1, stringsAsFactors=F, sep=",")[,(1:3)]
colnames(cols1) <- c("POP", "gen", "popSize")

avgRes <- Reduce("+", allres) / length(allres)


# dataframe with load values averaged acorss runs
d <- cbind(cols1, avgRes)
d$genbp <- max(d$gen) - d$gen
d <- d[d$genbp < 3e4 & d$genbp > 0,] # this discards first 30 k years 

firstgen <- max(d$genbp)
splitgen <- d$genbp[which(d$POP == "South")[1]]






#doPlot <- function() {

    #################################################
    ###### PANEL 1 PLOT POPULATION SIZES ############
    #################################################
    par(mar=c(0.5,4,1,4), oma=c(5,4,3.5,0))
    k <- d$POP == "North"

    ylims <- c(50, max(d$popSize) * 1.1)
    xlims <- c(max(d$genbp) + 1, 1)
    plot(y=d$popSize[k], x=d$genbp[k], log="xy", type="l", col=saola_colors["Northern"], ylim=ylims,xlim=xlims, 
    xaxt="n", cex.axis=cexaxis, ylab="Population size", cex.lab=cexlab, xpd=NA, xlab="", lwd=4)


    #k <- d$POP == "South"

    #segments(y0=sort(unique(d$popSize), decreasing=T)[2], y1=max(d$popSize[k]),  x0=d$genbp[which(d$POP == "South")[1]], lwd=4, col=saola_colors["Southern"])
    #lines(y=d$popSize[k], x=d$genbp[k] + 1, col=saola_colors["Southern"], lwd=4)

    #abline(v=splitgen, lty=2, lwd=3)

    #text(x=lastgen800, y=11000, labels="Population\nsplit",cex=1.8)

    #legend("topleft", col=saola_colors, legend=names(saola_colors), title="Population", lwd=4, bty="n", cex=2)

    arrows(x0=xlims[1], x1=xlims[2], code=2, y0=ylims[2] * 1.8, xpd=NA, lwd=3, length=0.1, angle=20)
    text(x=exp((log(xlims[2]) +  log(xlims[1])) / 2), y=ylims[2] * 2.75, xpd=NA, labels="Time", cex=2.5)
    ####################################################
    ###### PANEL 2 PLOT REALIZED GENETIC LOAD ##########
    ####################################################

    genloadcols <- 4:15 + 12
    d$realizedLoad <- rowSums(d[,genloadcols])
    ylims <- c(0, max(d$realizedLoad))
    xlims <- c(max(d$genbp) + 1, 1)

    # do empty plot
    plot(x=1, type="n", log="x",  col=saola_colors["Northern"],
         xlim=xlims, ylim=ylims,
         xaxt="n", cex.axis=cexaxis, ylab="Realized\ngenetic load", cex.lab=cexlab, xpd=NA, xlab="")

    lims <- c(0, 0.01, 0.05, 0.1, 0.25, 0.5)
    nbins <- length(lims)
    northpal <- rev(colorRampPalette(c(lighten(saola_colors["Northern"],0.2), darken(saola_colors["Northern"],0.6)))(nbins))
    #southpal <- rev(colorRampPalette(c(lighten(saola_colors["Southern"],0.2), darken(saola_colors["Southern"],0.6)))(nbins))

    lastgen <- 1
    plotgens <- d$genbp[d$genbp <= firstgen & d$genbp >= (lastgen)][k]
    htext <- c("[0, 0.01]", "(0.01, 0.05]", "(0.05, 0.1]", "(0.1, 0.25]", "(0.25, 0.5]", "(0.5, 1]")

    text(x=5500, y=mean(rowSums(d[d$genbp %in% plotgens,genloadcols])) + 0.1,
         labels="Dominance\ncoefficient", cex=1.4, font=2)

    for(i in 1:length(lims)) {
        
        cumsumload <- rowSums(d[,genloadcols[binStarts >= lims[i]]])
        k <- d$POP == "North"
        #lines(y=cumsumload[k], x=d$genbp[k] + 1, col=northpal[i], lwd=3, lty=1)

        # draw stacked polygon in generations previous to split
        polygon(x=c(firstgen, plotgens, lastgen), y=c(0, cumsumload[d$genbp %in% plotgens & k],0), col=northpal[i], lwd=0.5)
        midpoint <- mean(c(mean(cumsumload[d$genbp %in% plotgens]),
                           mean(ifelse(i<length(lims), rowSums(d[,genloadcols[binStarts >= lims[i + 1]]]), 0))))
        text(x=5500, y=midpoint,  htext[i], col=ifelse(i < 4, "white", "black"), cex=1.1)
        
    }

    #################################################
    ###### PANEL 3 PLOT PURGING GENETIC LOAD ##########
    #################################################

    genloadcols <- 4:15
    d$genLoad <- rowSums(d[,genloadcols])
    refLoad <- d$genLoad[1]
    ylims <- c(0, 1.2)
    xlims <- c(max(d$genbp) + 1, 1)#0000)

    # do empty plot
    plot(x=1, type="n", log="x",  col=saola_colors["Northern"],
         xlim=xlims, ylim=ylims,
         xaxt="n", cex.axis=cexaxis,
         ylab="",
        # ylab=expression("Purging\n(Gen. load"["t=ref", d$genbp[1]] ~ " / Gen. load"["t=x"] ~ ")"),
         cex.lab=cexlab, xpd=NA, xlab="", yaxt="n")
    axis(2, at=seq(0,1,0.2), cex.axis=cexaxis)
    title(ylab="Purging", cex.lab=cexlab, line=5.3,xpd=NA)
    title(ylab=expression("(Gen. load"["t=x", d$genbp[1]] ~ " / Gen. load"["t=ref"] ~ ")"), cex.lab=cexlab-0.5, line=3,xpd=NA)

#    text(x=5500, y=mean(rowSums(d[d$genbp %in% plotgens,genloadcols]))/refLoad + 0.12,
#         labels="Dominance\ncoefficient", cex=1.4)

    for(i in 1:length(lims)){
        
        cumsumload <- rowSums(d[,genloadcols[binStarts >= lims[i]]])
        k <- d$POP == "North"
        #lines(y=cumsumload[k]/refLoad, x=d$genbp[k] + 1, col=northpal[i], lwd=3, lty=1)

        
        # draw stacked polygon in generations previous to split
        polygon(x=c(firstgen, plotgens, lastgen), y=c(0, cumsumload[d$genbp %in% plotgens & k]/refLoad,0), col=northpal[i], lwd=0.5)
        midpoint <- mean(c(mean(cumsumload[d$genbp %in% plotgens]/refLoad),
                           mean(ifelse(i<length(lims), rowSums(d[,genloadcols[binStarts >= lims[i + 1]]])/refLoad, 0))))
        if(i < 3){
            text(x=5500, y=midpoint,   htext[i], col=ifelse(i < 4, "white", "black"), cex=1.1)
        }
        #else{
        # text(x=13000, y=rev(seq(-0.01, 0.25, length.out=4))[i-2],  htext[i], cex=1.1, xpd=NA, pos=4, xpd=NA)
        # arrows(x0=13500, x1= 8000, y0=rev(seq(-0.01, 0.25, length.out=4))[i-2], y1=midpoint, code=2, length=0.05, angle=20, xpd=NA)
        #}
        

        # draw stacked barplot with final load
        linek <- k & !(d$genbp %in% plotgens)
        lines(y=cumsumload[linek]/refLoad, x=d$genbp[linek] + 1, col=northpal[i], lwd=2, lty=1) }

    abline(v=d$genbp[1], lty=2, lwd=2)
    text(x=d$genbp[1]+100, y=1.1, labels="t=ref",cex=1.75, pos=4)

        # k <- d$POP == "South"
        # lines(y=cumsumload[k]/refLoad, x=d$genbp[k] + 1, col=southpal[i], lwd=3, lty=1)

        # draw stacked barplot with final load
        # segments(y0=rev(cumsumload[k]/refLoad)[1], x0=min(d$genbp + 1), x1= min(d$genbp + 1)-exp(-0.5), col=southpal[i], lty=3, lwd=2)
        # rect(ytop=rev(cumsumload[k]/refLoad)[1], ybottom=0, xleft=min(d$genbp + 1)-exp(-0.5), xright=min(d$genbp+1)- exp(0), col=southpal[i]) }

    #rect(ytop=1, ybottom=0.9,
    #     xleft=exp(c(0, 0.5, 1, 1.5, 2, 2.5)),
    #     xright=exp(c(0, 0.5, 1, 1.5, 2, 2.5) + 0.5),
    #     col=southpal, border="black")
    binlims <- c(lims, 1) # c(1, rev(lims))
    rect(ytop=1.01, ybottom=0.91,
         xright=exp(c(0, 0.5, 1, 1.5, 2, 2.5)),
         xleft=exp(c(0, 0.5, 1, 1.5, 2, 2.5) + 0.5),
         col=rev(northpal))
    text(x=exp(c(0, 0.5, 1, 1.5, 2, 2.5, 3)), y = 1.06, labels=rev(binlims), cex=1.15)
    text(x=exp(1.5), y=1.17, labels="Dominance coefficient bins", cex=1.3, font=2)


    arrows(x0=exp(3), x1=exp(0), y0=0.83, code=2, length=0.1, angle=20)
    text(x=exp(c(0.5, 2.5)), y=0.8, pos=1, 
         labels=rev(c("More\ndeleterious", "Less\ndeleterious")), cex=1.3)

    #text(x=exp(c(0.5, 2.5)), y=0.87, pos=1, 
    #     labels=c("Recessive\nmutations", "Additive\n(and dominant)\nmutations"), cex=1.3)


    axis(1, at=c(1,10,100,1000, 10000) + 1, labels=c(1,10,100,1000, 10000),  cex.axis=cexaxis)
    title(xlab="Generations before present (t)", xpd=NA, cex.lab=cexlab, xpd=NA)
    
#}




#doPlot()
