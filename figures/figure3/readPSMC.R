source("/home/lpq293/mypopgen/saola/analyses/psmc_het/scripts/readpsmc.R")


# read in genome-wide psmc data

files_allmuts <- c("/home/lpq293/mypopgen/saola/analyses/psmc_het/results/psmc_output/all_filters/2/9264_2.psmc",
                   "/home/lpq293/mypopgen/saola/analyses/psmc_het/results/psmc_output/all_filters/2/9176merged.psmc",
                   "/home/lpq293/mypopgen/saola/analyses/psmc_het/results/psmc_output/all_filters/2/dups1.psmc")


res_psmc <- lapply(files_allmuts, psmc.result,mu=mu,g=g)


names(res_psmc) <- gsub(".psmc", "", basename(files_allmuts))




# read in psmc in roh intersections and get cross coalescence



files_allmuts <- c(
    "/home/lpq293/mypopgen/saola/analyses/rohintersections_stuff/results_6dp_backup/psmc_output/Saola_PNgh_sites_nomultiallelics_noindels_6dp_3het_sample9176merged_rohs9264_2subtracted9176merged_fqtag.psmc",
    "/home/lpq293/mypopgen/saola/analyses/rohintersections_stuff/results_6dp_backup/psmc_output/Saola_PNgh_sites_nomultiallelics_noindels_6dp_3het_sample9176merged_rohsdups1subtracted9176merged_fqtag.psmc",
    "/home/lpq293/mypopgen/saola/analyses/rohintersections_stuff/results_6dp_backup/psmc_output/Saola_PNgh_sites_nomultiallelics_noindels_6dp_3het_sample9264_2_rohs9176mergedsubtracted9264_2_fqtag.psmc",
    "/home/lpq293/mypopgen/saola/analyses/rohintersections_stuff/results_6dp_backup/psmc_output/Saola_PNgh_sites_nomultiallelics_noindels_6dp_3het_sample9264_2_rohsdups1subtracted9264_2_fqtag.psmc",
    "/home/lpq293/mypopgen/saola/analyses/rohintersections_stuff/results_6dp_backup/psmc_output/Saola_PNgh_sites_nomultiallelics_noindels_6dp_3het_sampledups1_rohs9176mergedsubtracteddups1_fqtag.psmc",
    "/home/lpq293/mypopgen/saola/analyses/rohintersections_stuff/results_6dp_backup/psmc_output/Saola_PNgh_sites_nomultiallelics_noindels_6dp_3het_sampledups1_rohs9264_2subtracteddups1_fqtag.psmc"
)



res_allmuts <- lapply(files_allmuts, psmc.result,mu=mu,g=g)


names(res_allmuts) <- c("9176merged", "9176merged", "9264_2", "9264_2", "dups1", "dups1")



files_fakediploid <- c("/home/lpq293/mypopgen/saola/analyses/rohintersections_stuff/results_6dp/psmc_output/Saola_PNgh_sites_nomultiallelics_noindels_6dp_3het_pseudodiploid_fake_9264_2_9176merged.psmc",
                   "/home/lpq293/mypopgen/saola/analyses/rohintersections_stuff/results_6dp/psmc_output/Saola_PNgh_sites_nomultiallelics_noindels_6dp_3het_pseudodiploid_fake_dups1_9176merged.psmc",
                   "/home/lpq293/mypopgen/saola/analyses/rohintersections_stuff/results_6dp/psmc_output/Saola_PNgh_sites_nomultiallelics_noindels_6dp_3het_pseudodiploid_fake_9264_2_dups1.psmc")

res_fakediploid <- lapply(files_fakediploid, psmc.result,mu=mu,g=g)

names(res_fakediploid) <- c("9264_2-9176merged", "dups1-9176merged", "9264_2-dups1")





bwpop <- c("9264_2-9176merged", "dups1-9176merged")
whithpop <- "9264_2-dups1"

withpop <- 1 / rowMeans(sapply(res_allmuts, function(x) x$Ne))
#withpop <- 1 /  rowMeans(sapply(res_fakediploid[whithpop], function(x) x$Ne))
#withpop <- 1 /  rowMeans(sapply(list(res_allmuts, res_fakediploid[whithpop]), function(x) x$Ne))
betwpop <- 1 / rowMeans(sapply(res_fakediploid[bwpop], function(x) x$Ne))

withpop <- 1 / apply(sapply(res_allmuts, function(x) x$Ne), 1, mean)
betwpop <- 1 / apply(sapply(res_fakediploid[bwpop], function(x) x$Ne),1, mean)


crosscoalrate <- betwpop / withpop

alltimes <- sapply(append(res_allmuts, res_fakediploid[bwpop]), function(x) x$YearsAgo)
crosscoaltimes <- apply(alltimes,1, mean)



