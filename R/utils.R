calculate_sigmas <- function(sigma_min, sigma_max, n = 10) {
  sigma_min <- log10(sigma_min)
  sigma_max <- log10(sigma_max)

  10 ^ seq(sigma_min, sigma_max, length.out = n)
}

# define simpson integration function
simpson <- function(x, y) {
  if (length(x) != length(y)) {
    stop("x and y must be the same length")
  }
  if (length(x) < 3) {
    stop("x and y must have at least 3 elements")
  }
  if (any(diff(x) <= 0)) {
    stop("x must be strictly increasing")
  }
  if (any(abs(diff(diff(x))) > 1e-10)) {
    stop("x must be equally spaced")
  }

  h <- diff(x)
  len_y <- length(y)
  integral_value <- (h[1] / 3) * (y[1] + 4 * sum(y[seq(2, len_y, by = 2)]) + 2 * sum(y[seq(3, len_y - 1, by = 2)]) + y[len_y])
  integral_value
}


volume_square <- function(ch, points, sigma = 1) {
    n <- nrow(points)
    v <- c()
    sigma <- sigma / 2

    # Extract the points that determine the convex hull
    points_hull <- ch$p[ch$hull, ]
    points <- as.matrix(points)

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

        cube_inhull <- geometry::inhulln(ch, points_i)
        if (all(cube_inhull)) {
            vol <- (sigma * 2) ^ 3
        } else {
            convhull <- geometry::intersectn(points_hull, points_i)
            vol <- convhull$ch$vol
        }

        v <- c(v, vol)
    }
    v / n
}
