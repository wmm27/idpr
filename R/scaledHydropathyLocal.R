#' Calculate the Average Scaled Hydropathy of an Amino Acid Sequence
#'
#' This is used to calculate the scaled hydropathy of an amino acid
#'   sequence using a sliding window. The output is either a dataframe or graph
#'   showing the calculated scores for each window along the sequence.
#' @inheritParams sequenceCheck
#'
#' @param window a positive, odd integer. 7 by default.
#'   Sets the size of sliding window, must be an odd number.
#'   The window determines the number of residues to be analyzed and averaged
#'   for each position along the sequence.
#' @param plotResults logical value, TRUE by default.
#'   If \code{plotResults = TRUE} a plot will be the output.
#'   If \code{plotResults = FALSE} the output is a data frame with scores for
#'   each window analyzed.
#' @param proteinName character string with length = 1.
#'   optional setting to replace the name of the plot if hydropathy = TRUE.
#' @param ... any additional parameters, especially those for plotting.
#' @return see plotResults argument
#' @family scaled hydropathy functions
#' @seealso \code{\link{KDNorm}} for residue values.
#' @references Kyte, J., & Doolittle, R. F. (1982). A simple method for
#'   displaying the hydropathic character of a protein.
#'   Journal of molecular biology, 157(1), 105-132.
#' @export




scaledHydropathyLocal <- function(
  sequence,
  window = 9,
  plotResults = TRUE,
  proteinName = NA,
  ...) {

  seqVector <- sequenceCheck(
    sequence = sequence,
    method = "stop",
    outputType = "vector",
    supressOutputMessage = T)

  if ((window %% 2) == 0) {
    stop("Window must be an odd number")
  }

  if (!all(c(is.logical(plotResults)))) {
    stop("plotResults and centerResidue
         require logical values")
  }

  names(seqVector) <- NULL
  seqLength <- length(seqVector)
  numberResiduesAnalyzed <- seqLength - (window - 1)

  #--------
  positionVector <- ((window - 1) / 2 + 1): (seqLength - (window - 1) / 2)
  centerResidueVector <- seqVector[positionVector]
  windowVector <- rep(NA, numberResiduesAnalyzed)
  scoreVector <- rep(NA, numberResiduesAnalyzed)

  #----------
  #Analysis
  for (i in 1:numberResiduesAnalyzed) {

    windowBegining.i <- i
    windowEnd.i <- i + (window - 1)

    sequenceWindow <- seqVector[windowBegining.i:windowEnd.i]
    windowVector[i] <- paste0(sequenceWindow, collapse = "")

    windowValues <- KDNorm$V2[match(sequenceWindow, KDNorm$V1)]
    #match gets the positions of the letters in the kdnorm table.
    #then the values of the normalized table is gathered based on the positions,
    # and is stored into the windowValues var

    scoreVector[i] <- sum(windowValues) / window
    #This is the hydropathy for the residue based on the window size
  }

  windowDF <- data.frame(Position = positionVector,
                         Window = windowVector,
                         CenterResidue = centerResidueVector,
                         WindowHydropathy = scoreVector)
  #---------
  #Output
  if (plotResults) {

    if (!is.na(proteinName)) {
      plotTitle <- paste0("Measurement of Scaled Hydropathy in ", proteinName)
    } else {
      plotTitle <- "Measurement of Scaled Hydropathy"
    }
    meanScaledHydropathyValue <- meanScaledHydropathy(sequence = sequence,
                                                      roundScore = 3)
    plotSubtitle <- paste0("Window Size = ",
                           window,
                           " ; Average Scaled Hydropathy = ",
                           meanScaledHydropathyValue)
    gg <-  sequencePlot(
      position = windowDF$Position,
      property = windowDF$WindowHydropathy,
      hline = meanScaledHydropathyValue,
      dynamicColor = windowDF$WindowHydropathy,
      customColors = c("chocolate1", "skyblue3", "grey65"),
      customTitle = NA,
      propertyLimits = c(0, 1))

    gg <- gg + ggplot2::labs(title = plotTitle,
                             subtitle = plotSubtitle)
    return(gg)

  } else { #returns the DF
    return(windowDF)
  }
}
