#' Protein Scaled Hydropathy Calculations
#'
#' This is used to calculate the scaled hydropathy of an amino acid
#'   sequence for each residue in the sequence.
#'   Results are based on the
#'
#' @inheritParams sequenceCheck

#' @param plotResults logical value, FALSE by default.
#'   If \code{plotResults = TRUE} a plot will be the output.
#'   If \code{plotResults = FALSE} the output is a data frame for each residue.
#' @param proteinName character string with length = 1.
#'   optional setting to include the name in the plot title.
#' @param ... any additional parameters, especially those for plotting.
#' @return if \code{plotResults = TRUE}, a graphical representation data.
#'   Average is shown by the horizontal line.
#'   If \code{plotResults = FALSE}, a dataframe is reported
#'   with each amino acid and each residue value shown.
#'   Score for each residue shown in the column "Hydropathy".
#' @family scaled hydropathy functions
#' @seealso \code{\link{KDNorm}} for residue values.
#' @references Kyte, J., & Doolittle, R. F. (1982). A simple method for
#'   displaying the hydropathic character of a protein.
#'   Journal of molecular biology, 157(1), 105-132.
#' @export



scaledHydropathyGlobal <- function(
  sequence,
  plotResults = FALSE,
  proteinName = NA,
  ...) {

  seqCharacterVector <- sequenceCheck(
    sequence = sequence,
    sequenceName = proteinName,
    method = "Stop",
    outputType = "Vector",
    supressOutputMessage = T
  )

  if (!is.logical(plotResults)) {
    stop("plotResults must be a logical value")
  }

  seqLength <- length(seqCharacterVector)


  scoreVector <- KDNorm$V2[match(seqCharacterVector, KDNorm$V1)]

  hydropathyDF <- data.frame(Position = c(1:seqLength),
                             AA = seqCharacterVector,
                             Hydropathy = scoreVector)

  if (plotResults) {

    meanScaledHydropathy <- sum(hydropathyDF$Hydropathy) / seqLength
    meanScaledHydropathy <- round(meanScaledHydropathy, 3)


    if (!is.na(proteinName)) {
      plotTitle <- paste0("Scaled Hydropathy of ", proteinName)
    } else {
      plotTitle <- "Scaled Hydropathy"
    }
    plotSubtitle <- paste0("Average Scaled Hydropathy = ",
                           meanScaledHydropathy)

    gg <-  sequencePlot(
      position = hydropathyDF$Position,
      property = hydropathyDF$Hydropathy,
      hline = meanScaledHydropathy,
      dynamicColor = hydropathyDF$Hydropathy,
      customColors = c("chocolate1", "skyblue3", "grey65"),
      customTitle = NA,
      propertyLimits = c(0, 1))

    gg <- gg + ggplot2::labs(title = plotTitle,
                             subtitle = plotSubtitle)

    return(gg)
  } else {
    return(hydropathyDF)
  }
}
