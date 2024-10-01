#!/usr/bin/env Rscript

source("scripts/plotting.R")
source("scripts/sfs.R")

options(scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
sfs_path <- args[1]

sfs <- sfs_path %>%
    sfs_read_path %>%
    sfs_normalise

df <- tibble(
    tajimas_theta = sfs_tajimas_theta(sfs),
    wattersons_theta = sfs_wattersons_theta(sfs),
    wus_theta = sfs_wus_theta(sfs),
)

df %>%
    mutate(across(everything(), ~as.character(signif(., 5)))) %>%
    format_tsv %>%
    cat

#stop() # For debugging
sfs_path <- "results/giraffe/sfs/RothschildsGiraffe_Masai.sfs"
