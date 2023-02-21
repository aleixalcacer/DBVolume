get_samples <- function(points, n) {
    samples <- purrr::map(
        .x = 1:3,
        .f = function(i) {
            samples <- stats::runif(n = n, min = min(points[, i]), max = max(points[, i]))
            tibble::tibble(!!as.character(i) := samples)
        }
    ) %>% purrr::reduce(.f = ~ dplyr::bind_cols(.x, .y))
    stats::setNames(samples, c("x", "y", "z"))
}

f <- function(point, means, sd) {
    res <- 0
    p <- as.matrix(point)
    for (i in seq_len(nrow(means))) {
        m <- as.matrix(means[i, ])
        res <- res + prod(stats::dnorm(p, mean = m, sd = sd))
    }
    res / nrow(means)
    res
}

volume <- function(points, ch, means = means, sd = 0.1) {
    in_hull <- geometry::inhulln(ch, as.matrix(points))
    n <- length(in_hull)
    v <- c()
    for (i in 1:n) {
        if (in_hull[i]) {
            v <- c(v, f(points[i, ], means = means, sd = sd))
        } else {
            v <- c(v, 0)
        }
    }
    v / n
}

inhull <- function(p, ch) {
    geometry::inhulln(ch, p)
}

volume_AB <- function(points, chA, chB, chI, mA, mB, mAB, sd = 0.1) {
    ih1 <- inhull(as.matrix(points), chA)
    ih2 <- inhull(as.matrix(points), chB)
    in_hull <- mapply(function(x, y) {
        x || y
    }, ih1, ih2)

    n <- length(in_hull)

    v <- c()
    for (i in 1:n) {
        if (in_hull[i]) {
            v <- c(v, f(points[i, ], means = mAB, sd = sd))
        } else {
            v <- c(v, 0)
        }
    }
    v / n
}
