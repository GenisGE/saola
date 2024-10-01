suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(patchwork)
    library(purrr)
    library(readr)
    library(scales)
    library(tidyr)
})

base_theme <- theme_linedraw(base_size = 18) %+replace%
    theme(
        legend.justification = "right",
        legend.margin = margin(t = -5, r = 0, b = 0, l = 0, unit = "pt"),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid = element_line(colour = "grey80"),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        strip.background = element_blank(),
    )

theme_set(base_theme)
