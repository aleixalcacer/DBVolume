#'
#' @title Compute the optimal sigma for the volume of the intersection
#' @description Compute the optimal sigma for the volume of the intersection
#' @param pointsA A matrix of points
#' @param pointsB A matrix of points
#' @param epsilon The error between the estimated volume and the PAVO volume
#' to determine the optimal sigma
#' @param n The number of samples to use
#' @param n_reps The number of repetitions to use
#'
#' @return The optimal sigma for the volume of the intersection
#'
#' @importFrom rlang :=
#' @importFrom dplyr %>%
#' 
#' @export
optimal_sigma <- function(pointsA, pointsB, epsilon = 0.01, n = 1e4, n_reps = 5) {
    # Compute the CH of the intersection
    convhull <- geometry::intersectn(pointsA, pointsB)
    chA <- convhull$ch1
    chB <- convhull$ch2

    chI <- convhull$ch

    pointsU <- rbind(pointsA, pointsB)

    epsilon <- 0.01

    max_sigma <- 0
    s <- 0.01

    # Compute PAVO alpha-method
    pavo_info <- pavo::voloverlap(pointsA, pointsB)


    # Compute the maximum sigma
    while (TRUE || s < 100) {
        vI <- 0
        vU <- 0
        for (i in 1:n_reps) {
            samples <- get_samples(pointsU, n)

            sd <- s

            vI <- vI + sum(volume(samples, chI, means = pointsU, sd = sd)) / n_reps
            vU <- vU + sum(volume_AB(samples, chA, chB, chI, pointsA, pointsB, pointsU, sd = sd)) / n_reps
        }

        if (abs(vI / vU - pavo_info$vboth) < epsilon) {
            max_sigma <- s
            break
        }
        s <- s + 5
    }


    # Compute the value of sigma between 0 and the max sigma where the volume of the intersection is max.

    max_vI <- -1
    rand_sigma <- stats::runif(5, 0, max_sigma)
    for (s in rand_sigma) {
        vI <- 0
        vU <- 0

        n_reps <- 5
        for (i in 1:n_reps) {
            samples <- get_samples(pointsU, n)

            sd <- s

            vI <- vI + sum(volume(samples, chI, means = pointsU, sd = sd)) / n_reps
        }

        if (vI > max_vI) {
            max_sigma <- s
            max_vI <- vI
        }
    }

    return(s)
}
