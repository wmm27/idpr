#' Calculate the Mean Scaled Hydropathy
#'
#' This funciton utlilizes the scaledHydropathyGlobal() function and
#'   easily returns the averaged hydropathy as a numeric value.
#'
#' @inheritParams sequenceCheck
#' @param roundScore Number of decimals the score will be rounded to.
#'   NA by default.
#' @return A numeric value equal to the Mean Scaled Hydropathy.
#' @family scaled hydropathy functions
#' @seealso \code{\link{KDNorm}} for residue values.
#' @references Kyte, J., & Doolittle, R. F. (1982). A simple method for
#'   displaying the hydropathic character of a protein.
#'   Journal of molecular biology, 157(1), 105-132.
#' @export

meanScaledHydropathy <-
  function(sequence,
           roundScore = NA) {

    hydropDF <- scaledHydropathyGlobal(sequence = sequence,
                                       plotOutput = FALSE)
    seqLength <- nrow(hydropDF)
    totalHydrop <- sum(hydropDF$Hydropathy)
    avgHydrop <- totalHydrop / seqLength

    if (is.na(roundScore)) {
      return(avgHydrop)
    } else {
      avgHydrop <- round(avgHydrop, roundScore)
      return(avgHydrop)
    }
  }
