#!/usr/bin/env Rscript

source("scripts/sfs.R")
source("scripts/plotting.R")

args <- commandArgs(trailingOnly = TRUE)
sfs_path <- args[1]
samples_path <- args[2]
out_path <- args[3]

#stop()
## For debugging
#sfs_path <- "results/giraffe/sfs/RothschildsGiraffe_Reticulated.sfs"
#samples_path <- "results/giraffe/sfs/RothschildsGiraffe_Reticulated.sampleslist"
#out_path <- "test.pdf"

df <- sfs_path %>%
    sfs_read_path %>%
    sfs_fold %>%
    sfs_tidy %>%
    rename("count" = value) %>%
    mutate(
        freq = count / sum(count),
        segregating = !allele %in% c(0, max(allele)),
        show = freq <= max(freq * segregating),
        freq_bar = if_else(show, freq, NA_real_),
        freq_text = comma(freq, accuracy = 0.00001),
        is_max = freq * show == max(freq * show)
    )

title <- sfs_get_name(sfs_path)

subtitle <- paste("Folded SFS, total sites:", comma(sum(df$count), scale = 1e-6, suffix = "M"))

caption <- read_lines(samples_path, progress = FALSE) %>%
    split((seq_along(.) - 1) %/% 4) %>%
    lapply(paste, collapse = ", ") %>%
    unlist %>%
    paste(collapse = ",\n")


plot <- ggplot(df, aes(allele, freq_bar, label = freq_text)) +
    geom_col(
        colour = "black",
        alpha = 0.75
    ) +
    geom_text(
        aes(y = show * freq),
        angle = 90,
        size = 5,
        colour = ifelse(df$is_max, "white", "black"),
        hjust = ifelse(df$is_max, 1.2, -0.2),
    ) +
    scale_x_continuous(breaks = seq(0, max(df$allele), 2)) +
    scale_y_continuous(label = comma) +
    labs(
        x = "Alleles",
        y = "Frequency",
        title = title,
        subtitle = subtitle,
        caption = caption,
    ) +
    theme(
        axis.text.y = element_text(angle = 90, hjust = 0.5)
    )

ggsave(out_path, plot, width = 10, height = 8)
