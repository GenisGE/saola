# Rscript combineDepths.R group pathtochroms maxdepth
# args:
#      1 group group depth was done to (infiles will be assumed to be in relative folder depths/group/
#      2 pathtochrom path to file with list of chromosomes for which depth was estimated
#      3 maxdepth integer with -maxdepth used in angsd
#      4 outmain main folder within which output subfolders are located
library(data.table)


args <- commandArgs(trailingOnly=T)

group <- args[1]
chrsfile <- args[2]
maxdepth <- as.integer(args[3])
outmain <- args[4]

infolder <- paste0(outmain, "/depths/", group)
outfolder <- paste0(outmain, "/depths")

chrs <- scan(chrsfile, what="dl")

a <- numeric(maxdepth+1)

for(chr in chrs){

    l <- scan(paste0(infolder,"/",chr,"_", group,".depthGlobal"), what=3)
    a <- a + l
}

med <- median(rep.int(0:maxdepth, a))


outfile1 <- paste0(outfolder,"/all_", group,".depthGlobal")
write(a, outfile1, ncolumns=length(a), sep="\t")

outfile2 <- paste0(outfolder,"/", group,".depthMedian")
write(med, outfile2,ncolumns=1)
