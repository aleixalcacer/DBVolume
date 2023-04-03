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

volume_square <- function(ch, points, sigma = 1) {
    n <- nrow(points)
    v <- c()
    sigma <- sigma / 2
    
    # Extract the points that determine the convex hull
    points_hull <- ch$p[ch$hull, ]

    for (i in 1:n) {
        points_i <- matrix(0, nrow = 8, ncol = 3)

        points_i[1, ] <- points[i, ] + c(-sigma, -sigma, -sigma)
        points_i[2, ] <- points[i, ] + c(-sigma, -sigma, sigma)
        points_i[3, ] <- points[i, ] + c(-sigma, sigma, -sigma)
        points_i[4, ] <- points[i, ] + c(-sigma, sigma, sigma)
        points_i[5, ] <- points[i, ] + c(sigma, -sigma, -sigma)
        points_i[6, ] <- points[i, ] + c(sigma, -sigma, sigma)
        points_i[7, ] <- points[i, ] + c(sigma, sigma, -sigma)
        points_i[8, ] <- points[i, ] + c(sigma, sigma, sigma)

        convhull <- geometry::intersectn(points_hull, points_i)
        vol <- convhull$ch$vol
        v <- c(v, vol)
    }
    v / n
}
