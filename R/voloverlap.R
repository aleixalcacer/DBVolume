#'
#' @title Volume overlap estimation for two sets of points
#' @description Estimate the volume overlap between two sets of points
#' @param pointsA A matrix of points
#' @param pointsB A matrix of points
#' @param sigma The standard deviation of the Gaussian kernel
#'
#' @return A data frame with the volume overlap estimation
#'
#' @importFrom rlang :=
#' @importFrom dplyr %>%
#' 
#' @export
volume_overlap <- function(pointsA, pointsB, sigma = 1) {
    convhull <- geometry::intersectn(pointsA, pointsB)
    chA <- convhull$ch1
    chB <- convhull$ch2
    chI <- convhull$ch

    pointsU <- rbind(pointsA, pointsB)

    vI <- volume_square(chI, pointsU, sigma = sigma)
    vA <- volume_square(chA, pointsU, sigma = sigma)
    vB <- volume_square(chB, pointsU, sigma = sigma)

    vU <- vA + vB - vI

    vboth <- vI / vU
    vsmallest <- vI / min(vA, vB)

}


volume_overlap_normal <- function(pointsA, pointsB, sigma = 1, n = 1e4, n_reps = 5) {
    sd = sigma
    convhull <- geometry::intersectn(pointsA, pointsB)
    chA <- convhull$ch1
    chB <- convhull$ch2
    chI <- convhull$ch

    pointsU <- rbind(pointsA, pointsB)

    samples <- get_samples(pointsU, n)

    vI <- volume(samples, chI, means = pointsU, sd = sd)
    vU <- volume_AB(samples, chA, chB, chI, pointsA, pointsB, pointsU, sd = sd)

    # Params

    # Vsmall <- vI / min(v1, v2)

    Vboth <- sum(vI) / sum(vU)

    info <- data.frame(
        volCONVEX = pavo::voloverlap(pointsA, pointsB, type = "convex")$vboth,
        volALPHA = pavo::voloverlap(pointsA, pointsB, type = "alpha")$pboth,
        volDENSITY = Vboth
    )

    return(info)
}
