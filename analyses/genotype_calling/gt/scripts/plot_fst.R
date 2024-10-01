#!/usr/bin/env Rscript

source("scripts/plotting.R")
source("scripts/sfs.R")

get_names <- function(path) {
    path %>%
        basename %>%
        {strsplit(., "_")[[1]][2]} %>%
        {strsplit(., "\\.")[[1]][1]} %>%
        {strsplit(., "-")[[1]]}
}

read_fst <- function(path) {
    read_tsv(path, col_types = cols_only(hudsons_fst = col_double()), progress = FALSE) %>%
        pull(hudsons_fst)
}

args <- commandArgs(trailingOnly = TRUE)
plot_path <- args[1]
csv_path <- args[2]
stats <- args[-c(1, 2)]

df <- stats %>%
    map_dfr(
        ~tibble(
            population = list(get_names(.)),
            fst = read_fst(.)
        ) %>%
            unnest_wider(population, names_sep = "")
    ) %>%
    rename("i" = population1, "j" = population2) %>%
    mutate(fst_label = comma(fst, accuracy = 0.001))

plot <- ggplot(df, aes(i, j, fill = fst, label = fst_label)) +
    geom_tile() +
    geom_text() +
    scale_fill_distiller(
        palette = "Spectral",
        guide = "none",
    ) +
    labs(
        title = "Fst",
        subtitle = "Based on 2D SFS",
        x = NULL,
        y = NULL,
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(plot_path, plot, width = 12, height = 12)
write_csv(select(df, -fst_label), csv_path)

#stop() # Debug
#stats <- list.files("results/giraffe/sfs/", pattern = "^RothschildsGiraffe_.*-.*.stats", full.names = TRUE)
