

#source("/maps/projects/popgen/people/lpq293/saola/info/loadPopInfo.R")
#source("~/unicph/mypopgen/saola/info/loadPopInfo.R")
saola_colors <- c(saola_colors, "Admixed"="purple")
library(colorspace)

northpal <- colorRampPalette(c(lighten(saola_colors["Northern"],0.5), darken(saola_colors["Northern"],0.5)))(3)
southpal <- colorRampPalette(c(lighten(saola_colors["Southern"],0.5), darken(saola_colors["Southern"],0.5)))(3)
mixpal <- colorRampPalette(c(lighten(saola_colors["Admixed"],0.5), darken(saola_colors["Admixed"],0.5)))(3)

#cat("set up colors \n")
basedir <- "/home/lpq293/mypopgen/saola/analyses/simulations/nwf_sims/longer/perezmodel"
#basedir <- "~/unicph/mypopgen/saola/analyses/simulations/nwf_sims/"
scenarios <- c("north4", "north12", "north24",
                "south4", "south12", "south24",
                "admix4", "admix12", "admix24")

nruns <- 200



survival <- do.call("rbind", lapply(1:nruns, function(x) read.table(paste0(basedir, "/", x, "/survival_log.txt"), h=T)))

survivalProp <- colSums(survival == "Yes") / nrow(survival)


survivalRes <- data.frame(population = rep(c("Admixed", "Southern", "Northern"), each=3),
                            n_founders = rep(c(4, 12, 24), 3),
                             survivalprop = survivalProp,
                             extinctionprop = 1 - survivalProp,
                             nruns = nruns)

#cat("read in data\n")
#survivalRes <- data.frame(population = rep(c("Admixed", "Southern", "Northern"), each=3),
##                            n_founders = rep(c(4, 12, 24), 3),
   #                          prop = c(survivalProp, 1 - survivalProp),
  #                           survival = rep(c("yes", "no"), each=length(survivalProp)),
    #                         nruns = nruns)

ord <- c(sapply(c("Northern", "Southern", "Admixed"), function(x) which(survivalRes$population == x)))
survivalRes <- survivalRes[ord,]

cisSurvival <- apply(survivalRes[,c("survivalprop", "nruns")], 1, function(x) binom.test(x=x[1] * x[2], n=x[2], conf.level=0.95)$conf.int)

cisExtinction <- apply(survivalRes[,c("extinctionprop", "nruns")], 1, function(x) binom.test(x=x[1] * x[2], n=x[2], conf.level=0.95)$conf.int)


#f <- paste0(basedir, "/", x, "/saola_breeding_early.csv")

readLastEach <- function(f, scen = scenarios){

    d <- read.table(f, h=T, sep=",")
    return(data.frame(t(sapply(scen, function(x) d[which(d$POP==x)[sum(d$POP==x)],]))))

}

allres <- lapply(1:nruns, function(x) readLastEach(paste0(basedir, "/", x, "/saola_breeding_early.csv")))

 
boxres <- data.frame(scenario = rep(scenarios, each=200),
                    survival = c(as.matrix(survival[,scenarios])),
                    heterozygosity = c(do.call("rbind", lapply(allres, function(x) unlist(x$heterozygosity)))),
                    genetic_load = c(do.call("rbind", lapply(allres, function(x) unlist(x$genLoad)))),
                    realized_load = c(do.call("rbind", lapply(allres, function(x) unlist(x$realizedLoad))))
                    )


boxres$n_founders <- gsub("\\D", "", boxres$scenario)
boxres$pop <- gsub("[0-9]", "", boxres$scenario)
atpos <- c(c(1,2,4,5,7,8), c(1,2,4,5,7,8) + 10, c(1,2,4,5,7,8) + 20)

#set stuff for ploting
atpos <- c(c(1,2,3,4,5,6) + c(0.1, -0.1), c(1,2,3,4,5,6) + c(0.1, -0.1) + 7, c(1,2,3,4,5,6) + c(0.1, -0.1) + 14)
colpal <- c("#B3D6F7", "#2F546E", "#FFD9B2", "#7B5504", "#CC9FFF", "#560186")
colpal <- c("#2F546E","#B3D6F7",  "#7B5504", "#FFD9B2", "#560186", "#CC9FFF")
popord <- c("Northern", "Southern", "Admixed")
ord <- order(survivalRes$n_founders, sapply(survivalRes$population, function(x) which(popord == x)))


#cat("read to do second plot\n")


#doPlot2 <- function(){

par(mar=c(0.5,5.5,0.5,0.5))
xlims <- c(0, 20)
ylims <- c(0, 1.5)

#####################################
#### 1 PLOT SURVIVAL PROPORTIONS ####
#####################################
plot(1, type="n", xlab="", ylab="", xlim=xlims, ylim=ylims, xaxt="n", yaxt="n", bty="n")
rect(xleft = atpos-0.4, xright = atpos+0.4, 
        ybottom = 0, ytop = c(rbind(survivalRes$extinctionprop, survivalRes$survivalprop)[,ord]),
        col=colpal)
axis(2, at=seq(0,1,0.25), cex.axis=cexaxis)
title(ylab="Proportion\nof simulations", cex.lab=2.5, xpd=NA, adj=0, line=4)

segments(x0=atpos,
        y0 = c(rbind(cisExtinction[2,ord], cisSurvival[2,ord])),
        y1 = c(rbind(cisExtinction[1,ord], cisSurvival[1,ord])),
        lwd = 3)

arrows(x0=atpos,
        y0 = c(rbind(cisExtinction[2,ord], cisSurvival[2,ord])),
        y1 = c(rbind(cisExtinction[1,ord], cisSurvival[1,ord])),
        lwd = 3,
        code=3, angle=90, length=0.05)

legend(x=6, y=1.9, title="Simulation outcome", legend=c("Extinction", "Survival"), 
                fill=grey.colors(6)[c(1,6)], bty="n", cex=2.25, xpd=NA)

legend(x=13.5, y=1.9, title="Founding population", 
legend=c("Northern", "Southern", "50% each"), fill=saola_colors, 
bty="n", cex=2.25, xpd=NA)


################################
#### 2 PLOT REALIZED LOADS ####
################################
ylims <- range(boxres$realized_load) # c(0, max(boxres$realized_load))
a <- boxplot(boxres$realized_load ~ factor(boxres$survival, levels=c("No", "Yes")) + 
                                factor(boxres$pop, levels=c("north", "south", "admix")) + 
                                factor(boxres$n_founders, levels=c("4", "12", "24")),
        at = atpos, boxwex=0.75, col=colpal, outline=F, ann=F, xaxt="n", ylim=ylims, xlim=xlims, frame.plot=FALSE,
        cex.axis=cexaxis)

xidx <-  as.integer(interaction(factor(boxres$survival, levels=c("No", "Yes")),  
        factor(boxres$pop, levels=c("north", "south", "admix")),
        factor(boxres$n_founders, levels=c("4", "12", "24"))))
points(y=boxres$realized_load, x = jitter(atpos[xidx], 1), bg=rep(colpal, 3)[xidx], pch=21, cex=1.5)
title(ylab="Realized\ngenetic load", cex.lab=2.5, xpd=NA, line=3.5)

################################
#### 3 PLOT HETEROZYGOSITYS ####
################################
ylims <- c(0, max(boxres$heterozygosity))
a <- boxplot(boxres$heterozygosity ~ factor(boxres$survival, levels=c("No", "Yes")) + 
                                factor(boxres$pop, levels=c("north", "south", "admix")) + 
                                factor(boxres$n_founders, levels=c("4", "12", "24")),
        at = atpos, boxwex=0.75, col=colpal, outline=F, ann=F, xaxt="n",ylim=ylims, xlim=xlims, frame.plot=FALSE,
        cex.axis=cexaxis)

xidx <-  as.integer(interaction(factor(boxres$survival, levels=c("No", "Yes")),  
        factor(boxres$pop, levels=c("north", "south", "admix")),
        factor(boxres$n_founders, levels=c("4", "12", "24"))))
points(y=boxres$heterozygosity, x = jitter(atpos[xidx], 1), bg=rep(colpal, 3)[xidx], pch=21, cex=1.5)
title(ylab="Genome-wide\nheterozygosity",  cex.lab=2.5, xpd=NA, line=4)


#axis(1, at=c(1.5, 3.5,5.5) + rep(c(0, 8, 16), each=3), labels=rep(c("Northern", "Southern", "50% Northern\n50% Southern"), 3), las=2)
axis(1, tick=F, at=3 + c(0, 7, 14), labels=c("4 founders", "12 founders", "24 founders"), cex.axis=2, font=2)
#}

