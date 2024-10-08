

indir <- "/home/genis/saola/analyses/psmc_het/results6dp/psmc_output/all_filters/1"
outname <- ""
mu <- 1.2e-8
g <- 6


# read pscm output psmc.result() adapted from function in https://datadryad.org/stash/dataset/doi:10.5061/dryad.0618v

##-------Rescale the ith iteration result of PSMC, and make ready for plotting
# file: result file from PSMC
# i.iteration: the ith iteration
# mu: mutation rate
# s: bin size
# g: years per generation

psmc.result<-function(file,i.iteration=25,mu=1e-8,s=100,g=1)
{
	X<-scan(file=file,what="",sep="\n",quiet=TRUE)
	
	START<-grep("^RD",X)
	END<-grep("^//",X)
	
	X<-X[START[i.iteration+1]:END[i.iteration+1]]
	
	TR<-grep("^TR",X,value=TRUE)
	RS<-grep("^RS",X,value=TRUE)
	
	write(TR,"temp.psmc.result")
	theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
	N0<-theta0/4/mu/s 
	
	write(RS,"temp.psmc.result")
	a<-read.table("temp.psmc.result")
	Generation<-as.numeric(2*N0*a[,3]) # a[,3] is t_k
	Ne<-as.numeric(N0*a[,4]) #a[,4] is lambda_k
	
	file.remove("temp.psmc.result")
	
	n.points<-length(Ne)
	YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
		Generation[n.points])*g
	Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
		Ne[n.points])
	
	data.frame(YearsAgo,Ne)
}



files <-list.files(indir,".psmc",full.names=T)

res <- lapply(files, psmc.result,mu=mu,g=g)
names(res) <- gsub(".psmc", "",basename(files))


# THIS IS GENERALIZABLE, BUT WILL FAIL IF THERE ARE MORE THAN 9 INDIVIDUALS TO PLOT.
inds <- names(res)
cols <- RColorBrewer::brewer.pal("Set1", n=length(inds))

# THIS IS SPECIFIC FOR WARTHOGS!
#source("/home/genis/warthogs/info/loadPopInfo.R")
#inds <- names(res)
#pops <- c("Zimbabwe", "Namibia", "Zambia", "Ghana", "Tanzania", "Desert", "Ghana", "Ghana")
#names(pops) <- inds

#cols <- wartCols

res2 <- lapply(res, function(x) x[-((nrow(x)-8):nrow(x)),])





rm_last <- function(x){ ## FUNCITON TO CLEAN DATA, IS A BIT ARBITRARY MIGHT BE GOOD TO CHECK AGAIN SOME DAY?

    nes <- unique(x$Ne)
    rmv <- nes[(length(nes)-1):length(nes)]
    x[!x$Ne%in%rmv,]
}

res2 <- lapply(res,rm_last)

res2 <- res
bitmap(outname,width=12,height=6,res=300)
#ymax <- max(sapply(res2,function(x)max(x$Ne)))
ymax <- 5e4 # ARBITRARY THRESHOLD BECAUSE PLOTS USUALY LOOK NICE WITH THIS
par(mar=c(5,5,4,15)+0.1)
plot(type='l', x=res2[[inds[1]]]$YearsAgo, log='x', y=res2[[inds[1]]]$Ne, col = cols[1],lwd=3,
     xlab=sprintf("Years ago (mu=%.2e, g=%.1f)",mu,g),
     ylab="Effective population size",cex.lab=1.5, xaxt='n',
     xlim=c(0.5 * 10^4, 1.5 * 10^6), ylim=c(0, ymax), cex.axis=1.5)
#xpos <- seq(1.5 * 10^4, 1.5 * 10^6,by=10000)
xpos <- c(0.5e4,1e4, 1.5e4, 5e4, 1e5, 2e5, 4e5, 6e5,1e6, 1.5e6)
xpos <- c(seq(1e5, 1e6, by=1e5), seq(1e6, 1e7, by=1e6))
xlabs <- c(1e5, 5e5, 1e6, 5e6, 1e7)
xpos <- c(0.5e4,seq(1e4,1e5,by=1e4), seq(1e5, 1e6, by=1e5))
xlabs <- c(0.5e4,2e4, 1e5, 5e5,1e6)
axis(1, at=xpos, labels=F)
axis(1, at=xlabs, tick=F, labels= paste(as.integer(xlabs/1e3), "kya"), cex.axis=1.5)

for(i in 2:length(res2)) lines(x=res2[[inds[i]]]$YearsAgo, y=res2[[inds[i]]]$Ne,  col = cols[i],lwd=3)
legend(legend=inds,
       col=cols,
       lty=1,lwd=3,,border=NA, bty="n",cex=1.5, xpd=NA,
       y=ymax*0.8, x=1.75e6)
dev.off()



