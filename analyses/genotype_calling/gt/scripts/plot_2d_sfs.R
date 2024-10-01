#!/usr/bin/env Rscript

source("scripts/sfs.R")
source("scripts/plotting.R")

options(scipen = 999)

SFS_2D_COLOURS <- c(
    "#6e40aa", "#a83cb3", "#df40a1",
    "#ff507a", "#ff704e", "#f89b31",
    "#d2c934", "#aff05b", "#77f75c",
    "#63f07e", "#53d9ac", "#47b3d3"
)

args <- commandArgs(trailingOnly = TRUE)
sfs_path <- args[1]
out_path <- args[2]

df <- sfs_path %>%
    sfs_read_path %>%
    sfs_fold %>%
    sfs_tidy %>%
    drop_na %>%
    rename("count" = value) %>%
    mutate(
        freq = count / sum(count),
        segregating = !(i + j) %in% c(0, max(i) + max(j)),
        show = freq <= max(freq * segregating),
        freq_fill = if_else(show, freq, NA_real_),
        freq_text = ifelse(show, NA_character_, comma(freq, accuracy = 0.0001)),
        is_max = freq * show == max(freq * show)
    )

title <- sfs_get_name(sfs_path)

subtitle <- paste("Folded SFS, total sites:", comma(sum(df$count), scale = 1e-6, suffix = "M"))

names <- strsplit(title, "-")[[1]]

plot <- ggplot(df, aes(i, j, fill = freq_fill, label = freq_text)) +
    geom_tile() +
    geom_text() +
    scale_x_continuous(breaks = seq(0, max(df$i), 2)) +
    scale_y_continuous(breaks = seq(0, max(df$j), 2)) +
    scale_fill_stepsn(
        colours = SFS_2D_COLOURS,
        n.breaks = 10,
        nice.breaks = TRUE,
        trans = "log10",
        show.limits = TRUE,
        guide = guide_coloursteps(
            title.vjust = 0.85,
            frame.colour = "black",
            ticks.colour = "black",
        )
    ) +
    coord_cartesian(expand = FALSE) +
    labs(
        title = title,
        subtitle = subtitle,
        x = names[1],
        y = names[2],
        fill = NULL
    ) +
    theme(
        legend.key.width = unit(40, "mm"),
        legend.position = "bottom",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.margin = margin(t = -5, r = 40, b = 0, l = 0, unit = "pt"),
    )

ggsave(out_path, plot, width = 12, height = 12)

#stop()
## For debugging
#sfs_path <- "results/giraffe/sfs/RothschildsGiraffe_WestAfrican-SouthernAfrican.sfs"
#out_path <- "test.pdf"

