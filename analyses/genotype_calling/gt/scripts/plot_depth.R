#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
counts <- args[1]
plot_path <- args[2]

if (length(args) > 2) {
    width <- as.numeric(args[3])
    height <- as.numeric(args[4])
} else {
    width <- 10
    height <- 6
}

df <- read.table(
    counts,
    col.names = c("Depth", "Count"),
    colClasses = c("integer", "integer")
)

quantiles <- quantile(
    rep.int(df$Depth, df$Count),
    probs = c(0.25, 0.5, 0.75, 0.999)
)
cutoff <- quantiles["99.9%"]
quantiles <- quantiles[c("25%", "50%", "75%")]

quantile_df <- data.frame(
    Depth = unname(quantiles),
    Quantile = names(quantiles)
)

depth_breaks <- c(pretty(df$Depth[df$Depth < cutoff], 10), quantiles)
depth_breaks <- unname(sort(depth_breaks))

depth_plot <- ggplot(df, aes(Depth, Count)) +
    geom_col() +
    geom_vline(
        aes(xintercept = Depth),
        quantile_df,
        linetype = "dashed",
        colour = "orange"
    ) +
    geom_text(
        aes(x = Depth, y = Inf, label = Quantile),
        quantile_df,
        angle = 90,
        hjust = 1,
        vjust = 1.1
    ) +
    scale_x_continuous(limits = c(0, cutoff), breaks = depth_breaks) +
    scale_y_continuous(labels = scales::comma) +
    labs(
        title = "Depth distribution",
        subtitle = paste0("Sites with depth in top 0.1% not shown")
    ) +
    theme_minimal() %+replace%
    theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

ggsave(plot_path, depth_plot, width = width, height = height)
