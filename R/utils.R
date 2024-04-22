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
    points_hull <-  ch$p[unique(c(ch$hull)),]
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
        cube_hull <- geometry::convhulln(points_i, "n FA")
        hull_incubehull <- geometry::inhulln(cube_hull, points_hull)

        if (all(cube_inhull)) {
            vol <- (sigma * 2) ^ 3
        }
        # else if (all(!hull_incubehull)) {
        #     print(!hull_incubehull)
        #     cat("No intersection\n")
        #  vol <- 0
        # }
        else {
            convhull <- geometry::intersectn(points_hull, points_i)
            vol <- convhull$ch$vol
        }

        v <- c(v, vol)
    }
    v / n
}


plot_points <- function(points_01, points_02, pca=NULL) {

    points_u = rbind(points_01, points_02)

    if (is.null(pca)) {

        # Compute the PCA
        pca <- prcomp(points_u)
    }

    # Transform points_u using pca$rotation

    pca_2d <- as.data.frame((as.matrix(points_u) %*% pca$rotation)[,1:2])

    pca_2d$group <- c(rep("points_01", nrow(points_01)), rep("points_02", nrow(points_02)))

    # compute CH for each group
    pca_points_01 <- pca_2d[pca_2d$group == "points_01", ]
    pca_points_02 <- pca_2d[pca_2d$group == "points_02", ]


    ch_points_01 <- pca_points_01[chull(pca_points_01[,1:2]),]
    ch_points_02 <- pca_points_02[chull(pca_points_02[,1:2]),]

    # Plot the PCA in 2D

    p <- ggplot() +
        geom_point(data = pca_points_01, aes(x = PC1, y = PC2), color = "green") +
        geom_point(data = pca_points_02, aes(x = PC1, y = PC2), color = "red") +
        geom_polygon(data = ch_points_01, aes(x =PC1, y = PC2), fill = "green", alpha = 0.1) +
        geom_polygon(data = ch_points_02, aes(x = PC1, y = PC2), fill = "red", alpha = 0.1) +
        labs(x = "PC1", y = "PC2", color = "Group")


    p
}
