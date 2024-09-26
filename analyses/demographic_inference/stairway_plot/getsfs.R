

## sfs used /home/genis/saola/analyses/sfs_unfolded/results/notrans2/sfs/Northern4.sfs

## get number of sites and polymorphic sites
# in R

f <- "/home/genis/saola/analyses/sfs_unfolded/results/notrans2/sfs/Northern4.sfs"
sfs <- scan(f, what=.3, skip=1)
sum(sfs)
# 1179089655
cat(sfs[-c(1, length(sfs))], "\n")
# 51656.98 45871.93 39324.98 30601.06 29915.36 27864.05 25734.98 24460.81 23079.86 22536.33 21251.28 19709.54 19935.4 21783.37 21994.18



## generate folded sfs
# in R
fold<- function(sfs){
    mid <- ceiling(length(sfs)/2)

    folded <- sfs + c(rev(sfs)[1:(mid-1)], rep(0, mid))
    folded <- folded[1:mid]
    return(folded)
}

f <- "/home/genis/saola/analyses/sfs_unfolded/results/notrans2/sfs/Northern4.sfs"
sfs <- scan(f, what=.3, skip=1)
sfs <- fold(sfs)
sum(sfs)
# 1179089655
cat(sfs[-1], "\n")
# 73651.16 67655.3 59260.38 50310.6 51166.64 50400.39 48814.84 24460.81
