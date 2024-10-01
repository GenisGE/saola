suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

sfs_read_header <- function(str) {
    matches <- regmatches(str, gregexpr("[[:digit:]]+", str))
    as.numeric(unlist(matches))
}

sfs_read_flat <- function(str) {
    scan(text = str, quiet = TRUE)
}

sfs_read_path <- function(path) {
    lines <- readLines(path)

    stopifnot(length(lines) == 2)
    dims <- sfs_read_header(lines[1])
    flat <- sfs_read_flat(lines[2])

    n_dims <- length(dims)
    if (n_dims == 1) {
        stopifnot(dims == length(flat))
        sfs <- flat
    } else if (n_dims == 2) {
        stopifnot(prod(dims) == length(flat))
        sfs <- matrix(flat, nrow = dims[1], byrow = TRUE)
    } else {
        stop()
    }

    sfs
}

sfs_write <- function(sfs, file = stdout()) {
    dim <- paste(sfs_dims(sfs), collapse = "/")
    header <- paste("# Alleles:", dim)
    sfs <- paste(sfs_flatten(sfs), collapse = " ")
    cat(header, sfs, sep = "\n", file = file)
}

sfs_get_name <- function(path) {
    gsub(".sfs", "", strsplit(basename(path), "_")[[1]][2])
}

sfs_get_populations <- function(path) {
    strsplit(sfs_get_name(path), "-")[[1]]
}

sfs_dims <- function(sfs) {
    n_dims <- sfs_n_dims(sfs)
    if (n_dims == 1) {
        length(sfs)
    } else {
        dim(sfs)
    }
}

sfs_flatten <- function(sfs) {
    n_dims <- sfs_n_dims(sfs)
    if (n_dims == 1) {
        sfs
    } else {
        c(t(sfs))
    }
}

sfs_n_dims <- function(sfs) {
    stopifnot(is.atomic(sfs) & is.numeric(sfs))
    if (is.null(dim(sfs))) {
        1
    } else {
        length(dim(sfs))
    }
}

sfs_fold <- function(sfs) {
    n_dims <- sfs_n_dims(sfs)
    if (n_dims == 1) {
        n <- length(sfs)
        stopifnot(n %% 2 == 1)
        mid <- (n - 1) / 2 + 1
        left <- 1:(mid - 1)
        right <- (mid + 1):n
        c(sfs[left] + rev(sfs[right]), sfs[mid])
    } else {
        folded <- matrix(NA, nrow = nrow(sfs), ncol = ncol(sfs))
        n <- dim(sfs)
        mid <- sum(n) / 2
        for (i in seq_len(nrow(sfs))) {
            for (j in seq_len(ncol(sfs))) {
                if (i + j < mid || (i + j == mid && i < j)) {
                    folded[i, j] <- sfs[i, j] + sfs[n[1] - i + 1, n[2] - j + 1]
                }
            }
        }
        folded
    }
}

sfs_normalise <- function(sfs, na.rm = FALSE) {
    sfs / sum(sfs, na.rm = na.rm)
}

sfs_is_normalised <- function(sfs) {
    isTRUE(all.equal(sum(sfs), 1))
}

sfs_marginalise <- function(sfs) {
    stopifnot(sfs_n_dims(sfs) == 2)
    list(rowSums(sfs), colSums(sfs))
}

sfs_tajimas_theta <- function(sfs) {
    stopifnot(sfs_n_dims(sfs) == 1)

    sfs <- sfs_normalise(sfs)
    n <- length(sfs) - 1

    seg_i <- 1:(n - 1)
    seg <- sfs[seg_i + 1]

    sum(seg_i * (n - seg_i) * seg) / choose(n, 2)
}

sfs_wus_theta <- function(sfs) {
    stopifnot(sfs_n_dims(sfs) == 1)

    sfs <- sfs_normalise(sfs)
    n <- length(sfs) - 1

    seg_i <- 1:(n - 1)
    seg <- sfs[seg_i + 1]

    sum(seg_i^2 * seg) / choose(n, 2)
}

sfs_wattersons_theta <- function(sfs) {
    stopifnot(sfs_n_dims(sfs) == 1)

    sfs <- sfs_normalise(sfs)
    n <- length(sfs) - 1

    seg_i <- 1:(n - 1)
    seg <- sfs[seg_i + 1]

    a <- sum(1 / seg_i)
    sum(seg) / a
}

sfs_f2 <- function(sfs) {
    stopifnot(sfs_n_dims(sfs) == 2)
    sfs <- sfs_normalise(sfs)

    n1 <- nrow(sfs) - 1
    n2 <- ncol(sfs) - 1
    p1 <- 0:n1 / n1
    p2 <- 0:n2 / n2

    f <- 0
    for (i in seq_along(p1)) {
        for (j in seq_along(p2)) {
            add <- (p1[i] - p2[j])^2 * sfs[i, j]
            f <- f + add
        }
    }

    f
}

sfs_hudsons_fst <- function(sfs) {
    stopifnot(sfs_n_dims(sfs) == 2)

    mat <- sfs
    mat[c(1, length(mat))] <- 0
    mat <- sfs_normalise(mat)

    n1 <- nrow(mat) - 1
    n2 <- ncol(mat) - 1
    p1 <- 0:n1 / n1
    p2 <- 0:n2 / n2

    f_num <- function(p1, p2) {
        q1 <- 1 - p1
        q2 <- 1 - p2
        (p1 - p2)^2 - p1 * q1 / (n1 - 1) - p2 * q2 / (n2 - 1)
    }
    num <- outer(p1, p2, f_num)

    f_denom <- function(p1, p2) {
        q1 <- 1 - p1
        q2 <- 1 - p2
        p1 * q2 + p2 * q1
    }
    denom <- outer(p1, p2, f_denom)

    sum(mat * num) / sum(mat * denom)
}

sfs_tidy <- function(sfs) {
    n_dims <- sfs_n_dims(sfs)
    if (n_dims == 1) {
        sfs_tidy_1d(sfs)
    } else {
        sfs_tidy_2d(sfs)
    }
}

sfs_tidy_1d <- function(sfs) {
    tibble(allele = seq_along(sfs) - 1, value = sfs)
}

sfs_tidy_2d <- function(sfs) {
    expand_grid(j = seq_len(ncol(sfs)) - 1, i = seq_len(nrow(sfs)) - 1) %>%
        mutate(value = c(sfs))
}
