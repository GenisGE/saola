source("/home/lpq293/mypopgen/saola/info/loadPopInfo.R")
library(UpSetR)

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


classify_depleted <- function(d, mindiv=0){
    
    res <- table(apply(d[,c("avg_pi.northern", "avg_pi.central", "avg_pi")] <= mindiv, 1, function(x) paste(as.integer(x), collapse=" ")))# / nrow(d)
    return(res)
}


get_res <- function(inprefix){

    categories <- c("0 0 0", "0 1 0", "1 0 0", "1 1 0", "1 1 1")
    d <- read_merge_tables(inprefix)
    res <- classify_depleted(d)

    categories <- c("0 0 0", "0 1 0", "1 0 0", "1 1 0", "1 1 1")

    # if a potential category is missing, set to 0
    if(any(! categories %in% names(res))) res[categories[! categories %in% names(res)]] <- 0
    res <- res[categories]

    # force order to be consistent
    return(res)

}


inprefix <- "/home/lpq293/mypopgen/saola/analyses/theta_win/results_matchedsamples/pixy_out/output_maxmiss0.5_win100000"

d <- read_merge_tables(inprefix)
res <- classify_depleted(d)
res2 <- res[-1]
names(res2) <- c("Southern", "Northern", "Southern&Northern", "Southern&Northern&Combined")

#res <- get_res(inprefix)
outpng <- "upset_depleteddivwin_100kb.png"
#bitmap(outpng, w=6,h=6,res=300)
png(outpng, w=6, h=6, unit="in", res=300)
upset(fromExpression(res2),
      main.bar.color=c(saola_colors[2:1], "darkblue", "darkred"),
      sets.bar.color=c(saola_colors[2:1], "black"),
      mainbar.y.label="Number of 100kb windows\nwithout genetic diversity",
      sets.x.label="Windows without\ngenetic diversity",
      text.scale=c(2,2,1.5,1.1,1.8,2.2) * 1.2,set_size.show=F)
dev.off()



outpdf <- "upset_depleteddivwin_100kb.pdf"
#bitmap(outpng, w=6,h=6,res=300)
pdf(outpdf, w=6, h=6)#, unit="in", res=300)
par(oma=c(1,2,0,0), xpd=NA)
upset(fromExpression(res2),
      main.bar.color=c(saola_colors[2:1], "darkblue", "darkred"),
      sets.bar.color=c(saola_colors[2:1], "black"),
      mainbar.y.label="Number of 100kb windows\nwithout genetic diversity",
      sets.x.label="Windows without\ngenetic diversity",
      text.scale=c(2,2,1.5,1.1,1.8,2.2) * 1.2,set_size.show=F)
dev.off()

