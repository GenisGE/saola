

### IMPLEMENT ESTIMATOR OF VARIANCE WITH UNEQUAL JACKKNIFE FROM busing et al 1999 Statisticn and computing

sfspre <- "../do_sfs/results/2dsfs/Southern-Northern"

readSumGlobal <- function()
    sum(scan(paste0(sfspre, ".sfs"), skip=1))

readSumSplits <- function(){

    a <- readLines((paste0(sfspre, "_split.sfs")))[seq(2, 100,2)]
    res <- sapply(a, function(x) sum(as.numeric(unlist(strsplit(x, " ")))))
    names(res) <- NULL
    return(res)
}


get_jk_est <- function(est_global, est_splits, n, ms, g=50)
    g * est_global - sum((1 - ms / n) * est_splits)


get_pseudo_est <- function(est_global, est_splits, n, ms, g=50){
    hs <- n / ms
    hs * est_global - (hs -1) * est_splits
}

get_jk_sd <- function(est_global, est_splits, n, ms, g=50){

    hs <- n / ms
    jk_est <- get_jk_est(est_global, est_splits, n, ms, g)
    jk_pseudo <- get_pseudo_est(est_global, est_splits, n, ms, g)
    jk_var <- (1/g) * sum((1/(hs - 1)) * ((jk_pseudo - jk_est) ^ 2))
    jk_se <- sqrt(jk_var)
    return(jk_se)
}
 
# number of sites in sfs
n <- readSumGlobal()

# number of sites removes for each split from the sfs
m_js <- readSumSplits()


maxlike <- read.table("../runfsc/model1/max_like_run/38/38.bestlhoods", h=T)
jk_ests <- read.table("model1_jacknife50.params", h=T)

sds <- sapply(1:(ncol(jk_ests) - 2),
              function(x) get_jk_sd(est_global=maxlike[,x], est_splits=jk_ests[,x], n=n, ms=m_js, g=50)) 
#sds2 <- sapply(1:(ncol(jk_ests) - 2),
#              function(x) get_jk_sd2(est_global=maxlike[,x], est_splits=jk_ests[,x], n=n, ms=m_js, g=50)) 
 
res <- t(rbind(maxlike[,1:(ncol(jk_ests) - 2)], sds))

outtsv <- "model1_jackknife_se.tsv"
write.table(res, outtsv, col.names=c("maxres", "se"), quote=F, row.names=T, sep="\t")

if(FALSE){

    # double check implementation using the second part of sd estiamtor in Busing et al  1999
    est_global <- maxlike[,i]
    est_splits <- jk_ests[,i]


    get_jk_sd2 <- function(est_global, est_splits, n, ms, g=50){

        hs <- n / ms
        jk_est <- get_jk_est(est_global, est_splits, n, ms, g)
        jk_var <- (1 / g) * sum((1 / (hs-1)) * (hs * est_global - (hs - 1) * est_splits - g * est_global + sum((1-(ms/n)) * est_splits)) ^ 2)
        #jk_var <- (1/g) * sum((1/(hs - 1)) * ((est_splits - jk_est) ^ 2))
        jk_se <- sqrt(jk_var)
        return(jk_se)
    }

    sds2 <- sapply(1:(ncol(jk_ests) - 2),
                function(x) get_jk_sd2(est_global=maxlike[,x], est_splits=jk_ests[,x], n=n, ms=m_js, g=50)) 
}