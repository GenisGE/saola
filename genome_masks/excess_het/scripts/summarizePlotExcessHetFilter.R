# Rscript summarizePlotExcessHetFilter.R inprefix outprefix badbed allbed minF minLRT
# inprefix prefix of pcangsd output files
# prefix for two output files (table with info per scaffold and big pdf with many plots
# allbed bedfile with whole chromosme (so 0 to chrsize)
# badbed bed with selected bad coordinates
# minF minimum F to consider site excess heterozygoisty
# minLRT minimum lrt test for significnat excess heterozygosity

library(RcppCNPy)
library(scales)

args <- commandArgs(trailingOnly=T)

inprefix <- args[1]
outprefix <- args[2]
badbedf <- args[3]
allbedf <- args[4]
minF <- args[5]
minLRT <- args[6]


badbed <- read.table(badbedf, col.names=c("chr", "start", "end"), stringsAsFactors=F)
allbed <- read.table(allbedf, col.names=c("chr", "start", "end"), stringsAsFactors=F)


estF <- npyLoad(paste0(inprefix,".inbreed.sites.npy"))
site <- scan(paste0(inprefix,".sites"),what="df")
lrt <- npyLoad(paste0(inprefix,".lrt.sites.npy"))
chr <- sub("_[0123456789]+$","",x=site)
pos <- as.integer(sub(".*_","",x=site))


meanF <- tapply(estF,chr,mean)
lenF <- tapply(estF,chr,length)

#maxPos <- tapply(pos,chr,max)


#lost <- sum(as.numeric(badbed$end-badbed$start))/sum(as.numeric(allbed$end))


#names(allbed) <- names(badbed)

# DO SUMMARY TABLE INDICATING SOME SUMMARY FILTER RESULTS PER SCAFFOLDS, MIGHT BE USED TO FILTER OUT WHOLE SCAFFOLDS
lenghtbad <- rep(0, nrow(allbed)) # this make sure chromosomes that have 0 bad regions are still represented as 0
lenghtbad[allbed$chr %in% unique(badbed$chr)] <- tapply(badbed$end-badbed$start, badbed$chr, sum)

Nbadreg <- sapply(allbed$chr, function(x) sum(badbed$chr==x))


summarydf <- data.frame(chr=allbed$chr, length=allbed$end,meanF=meanF,
                        proportionBad=lenghtbad/allbed$end,
                        Nbadregions=Nbadreg)

outtable <- paste0(outprefix, "_summary.tsv")

write.table(summarydf, outtable,
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


### DO BIG PDF WITH LOCAL VALUES AND FILTER
splitChrs <- function(c, s=10e6, bed=allbed){

    cat("doing chr",c,"\n")
    i <- which(bed[,1]==c)
    
    if(bed[i,3] < s) return(rep(paste(c,1,sep="_"), length(pos[chr==c])))
                            
    idx <- cut(x=pos[chr==c], breaks=seq(allbed[i,2], allbed[i,3], length.out=ceiling(allbed[i,3]/s)),
               labels=FALSE)

    

    return(paste(c, idx,sep="_"))

}
# split chromosomes in 10 mb fragments for better visualization
splitt <- unlist(sapply(allbed$chr, splitChrs))

meanF2 <- tapply(estF, splitt, mean)

w <- names(sort(meanF2))
k <- c(1:5, sort(sample(6:length(w), 15))) # keep to plot only 5 worst (most mean F scaffolds/segments plus 15 random ones)
w <- w[k]

titles <- paste(w,c(paste(1:5, "most negative mean F value region"), paste("random region", 1:15)))

outbigpdf <- paste0(outprefix, "_plotsBig.pdf")

# make huge pdf with all scaffolds
pdf(outbigpdf)
for(i in 1:length(w)){

    scaf <- w[i]
    ttitle <- titles[i]

    p<-pos[splitt==scaf]
    F<-estF[splitt==scaf]
    l <- lrt[splitt==scaf]

    plot(p,F,pch=4,col=ifelse(F< 1 & l>24, "goldenrod","grey"),lwd=2,ylab="F",xlab="position", main=ttitle,cex.lab=1.8,cex.axis=1.8,cex.main=1.8)

    win <- 200
    fun<-function(x,w) ( cumsum(as.numeric(x))[-c(1:w)]-rev(rev(cumsum(as.numeric(x)))[-c(1:w)])) / w
    lines(fun(p,win),fun(F,win))


    c <- sub("_[0123456789]+$","",x=scaf)

    
    start <- badbed[badbed$chr==c,]$start
    end <- badbed[badbed$chr==c,]$end
    abline(v=start, col="darkred",lwd=1)
    abline(v=end, col="darkred",lwd=1)

    #segments(y0=-1, y1=1, x0=start, x1=end, col="darkred")
    #segments(y0=1, y1=-1, x0=start, x1=end, col="darkred")

    rect(xleft=start,xright=end,ybottom=-1.5,ytop=1.5,col=scales::alpha("grey",0.5), border =NA)
}
dev.off()

