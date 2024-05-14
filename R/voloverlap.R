#'
#' @title Volume overlap estimation for two sets of points
#' @description Estimate the volume overlap between two sets of points
#' @param points_a A matrix of points
#' @param points_b A matrix of points
#' @param sigma The standard deviation of the Gaussian kernel. If NULL, the
#'  volume overlap is estimated for a range of sigmas.
#' @param n_sigma The number of sigmas to use for estimating the volume overlap.
#' Only used if sigma is NULL.
#'
#' @return A data frame with the volume overlap estimation
#'
#' @importFrom rlang :=
#' @importFrom dplyr %>%
#'
#' @export
volume_overlap <- function(points_a, points_b, sigma = NULL, n_sigma = 5) {
    convhull <- geometry::intersectn(points_a, points_b)
    ch_a <- convhull$ch1
    ch_b <- convhull$ch2
    ch_i <- convhull$ch
    points_u <- rbind(points_a, points_b)

    if (is.null(sigma)) {
        sigma_max <- (max(points_u) - min(points_u)) * 2
        sigma_min <- 1e-4
        sigmas <- calculate_sigmas(sigma_min, sigma_max, n = n_sigma)
    } else if (is.numeric(sigma)) {
      sigmas <- c(sigma)
    } else if (is.vector(sigma)) {
      sigmas <- sigma
    } else {
      stop("sigma must be a numeric or a vector")
    }

    vboth <- c()
    for (sigma in sigmas) {
        if (ch_i$vol == 0) {
            v_i_i <- 0
        } else {
            v_i_i <- sum(volume_square(ch_i, points_u, sigma = sigma))
        }
        v_a_i <- sum(volume_square(ch_a, points_u, sigma = sigma))
        v_b_i <- sum(volume_square(ch_b, points_u, sigma = sigma))
        v_u_i <- v_a_i + v_b_i - v_i_i
        vboth_i <- v_i_i / v_u_i
        # vsmallest_i <- v_i_i / min(v_a_i, v_b_i)

        # cat("v_i_i: ", v_i_i, "\n")
        # cat("v_a_i: ", v_a_i, "\n")
        # cat("v_b_i: ", v_b_i, "\n")


        vboth <- c(vboth, vboth_i)
    }

    vboth_sd <- sd(vboth)
    vboth_mean <- mean(vboth)

    df <- list(
        vboth = vboth_mean,
        vboth_sd = vboth_sd,
        vboth_cv = vboth_sd / vboth_mean,
        vboth_extended = vboth,
        sigmas = sigmas
    )

    return(df)
}
