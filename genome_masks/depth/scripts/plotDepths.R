# Rscript plotDepths.R outmain group1 group2 (...)
# jsut takes main folder and then id of groups that depths was done with, assumes there exists dist and median files for those groups

args <- commandArgs(trailingOnly=T)

outmain <- args[1]
groups <- args[2:length(args)]

ngroups <- length(groups)

infolder <- paste0(outmain, "/depths")

outfile <- paste0(outmain, "/plots/depthFilter.png")

bitmap(outfile, width=6, height=3*ngroups, res=300)
par(mfrow=c(ngroups,1))
for(group in groups){

    
    a <- scan(paste0(infolder,"/all_", group, ".depthGlobal"), what=1)
    med <- scan(paste0(infolder,"/", group, ".depthMedian"), what=1)

    low_thres <- med * 0.5
    high_thres <- med * 1.5
    
    max <- med * 2
    
    a <- c(a[1:max], sum(a[(max+1):(length(a))]))

    barplot(a, names.arg=0:max, main=paste("Depth distribution", group), space=0, border=NA)
    abline(v=med, col=1, lty=1)
    abline(v=low_thres, col=2, lty=2)
    abline(v=high_thres, col=2, lty=2)
}
dev.off()
