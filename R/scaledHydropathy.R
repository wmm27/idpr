#' Calculate the Average Scaled Hydropathy of an Amino Acid Sequence
#'
#' This is used to calculate the scaled hydropathy of an amino acid
#'   sequence using a sliding window. The output is either a data frame or graph
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
#' @examples
#' #Amino acid sequences can be character strings
#' aaString <- "ACDEFGHIKLMNPQRSTVWY"
#' #Amino acid sequences can also be character vectors
#' aaVector <- c("A", "C", "D", "E", "F",
#'               "G", "H", "I", "K", "L",
#'               "M", "N", "P", "Q", "R",
#'               "S", "T", "V", "W", "Y")
#' #Alternativly, .fasta files can also be used by providing
#'   ##The path to the file as a character string.
#'
#' exampleDF <- scaledHydropathyLocal(aaString,
#'                                    plotResults = FALSE)
#' head(exampleDF)
#'
#' exampleDF <- scaledHydropathyLocal(aaVector,
#'                                    plotResults = FALSE)
#' head(exampleDF)
#'
#' #Changing window will alter the number of residues analyzed
#' exampleDF_window3 <- scaledHydropathyLocal(aaString,
#'                                            window = 3,
#'                                            plotResults = FALSE)
#' head(exampleDF_window3)
#' exampleDF_window15 <- scaledHydropathyLocal(aaString,
#'                                             window = 15,
#'                                             plotResults = FALSE)
#' head(exampleDF_window15)
#'
#' #plotResults = TRUE will output a ggplot
#' \dontrun{
#'   scaledHydropathyLocal(aaString,
#'                         plot = T)
#'
#' #since it is a ggplot, you can change or annotate the plot
#'  gg <- scaledHydropathyLocal(aaVector,
#'                              window = 3,
#'                              plot = T)
#'  gg <- gg + ggplot2::ylab("Local Hydropathy")
#'   gg <- gg + ggplot2::geom_text(data = exampleDF_window3,
#'                                ggplot2::aes(label = CenterResidue,
#'                                             y = WindowHydropathy + 0.1))
#'  plot(gg)
#' }




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
    supressOutputMessage = TRUE)

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
  for (i in seq_len(numberResiduesAnalyzed)) {

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



#' Protein Scaled Hydropathy Calculations
#'
#' This is used to calculate the scaled hydropathy of an amino acid
#'   sequence for each residue in the sequence.
#'   The output is either a dataframe or graph
#'   showing the matched scores for each residue along the sequence.
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
#' @examples
#' #Amino acid sequences can be character strings
#' aaString <- "ACDEFGHIKLMNPQRSTVWY"
#' #Amino acid sequences can also be character vectors
#' aaVector <- c("A", "C", "D", "E", "F",
#'               "G", "H", "I", "K", "L",
#'            "M", "N", "P", "Q", "R",
#'            "S", "T", "V", "W", "Y")
#' #Alternativly, .fasta files can also be used by providing
#' ##The path to the file as a character string
#'
#' exampleDF <- scaledHydropathyGlobal(aaString,
#'                                     plotResults = FALSE)
#' head(exampleDF)
#'
#' exampleDF <- scaledHydropathyGlobal(aaVector,
#'                                     plotResults = FALSE)
#' head(exampleDF)
#'
#' #plotResults = TRUE will output a ggplot
#' \dontrun{
#'   scaledHydropathyGlobal(aaString,
#'                          plot = T)
#'
#'   #since it is a ggplot, you can change or annotate the plot
#'   gg <- scaledHydropathyGlobal(aaVector,
#'                                plot = T)
#'   gg <- gg + ggplot2::ylab("Local Hydropathy")
#'   gg <- gg + ggplot2::geom_text(data = exampleDF,
#'                                 ggplot2::aes(label = AA,
#'                                              y = Hydropathy + 0.1))
#'   plot(gg)
#' }

scaledHydropathyGlobal <- function(
  sequence,
  plotResults = FALSE,
  proteinName = NA,
  ...) {

  seqCharacterVector <- sequenceCheck(
    sequence = sequence,
    method = "stop",
    outputType = "vector",
    supressOutputMessage = TRUE
  )

  if (!is.logical(plotResults)) {
    stop("plotResults must be a logical value")
  }

  seqLength <- length(seqCharacterVector)


  scoreVector <- KDNorm$V2[match(seqCharacterVector, KDNorm$V1)]

  hydropathyDF <- data.frame(Position = seq_len(seqLength),
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

#' Protein Scaled Hydropathy Calculations
#'
#' This is used to calculate the scaled hydropathy of an amino acid
#'   sequence for each residue in the sequence.
#'   The output is either a data frame or graph
#'   showing the matched scores for each residue along the sequence.
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
#'   If \code{plotResults = FALSE}, a data frame is reported
#'   with each amino acid and each residue value shown.
#'   Score for each residue shown in the column "Hydropathy".
#' @family scaled hydropathy functions
#' @seealso \code{\link{KDNorm}} for residue values.
#' @references Kyte, J., & Doolittle, R. F. (1982). A simple method for
#'   displaying the hydropathic character of a protein.
#'   Journal of molecular biology, 157(1), 105-132.
#' @export
#' @examples
#' #Amino acid sequences can be character strings
#' aaString <- "ACDEFGHIKLMNPQRSTVWY"
#' #Amino acid sequences can also be character vectors
#' aaVector <- c("A", "C", "D", "E", "F",
#'               "G", "H", "I", "K", "L",
#'            "M", "N", "P", "Q", "R",
#'            "S", "T", "V", "W", "Y")
#' #Alternativly, .fasta files can also be used by providing
#' ##The path to the file as a character string
#'
#' exampleDF <- scaledHydropathyGlobal(aaString,
#'                                     plotResults = FALSE)
#' head(exampleDF)
#'
#' exampleDF <- scaledHydropathyGlobal(aaVector,
#'                                     plotResults = FALSE)
#' head(exampleDF)
#'
#' #plotResults = TRUE will output a ggplot
#' \dontrun{
#'   scaledHydropathyGlobal(aaString,
#'                          plot = T)
#'
#'   #since it is a ggplot, you can change or annotate the plot
#'   gg <- scaledHydropathyGlobal(aaVector,
#'                                plot = T)
#'   gg <- gg + ggplot2::ylab("Local Hydropathy")
#'   gg <- gg + ggplot2::geom_text(data = exampleDF,
#'                                 ggplot2::aes(label = AA,
#'                                              y = Hydropathy + 0.1))
#'   plot(gg)
#' }

scaledHydropathyGlobal <- function(
  sequence,
  plotResults = FALSE,
  proteinName = NA,
  ...) {

  seqCharacterVector <- sequenceCheck(
    sequence = sequence,
    method = "stop",
    outputType = "vector",
    supressOutputMessage = TRUE)

  if (!is.logical(plotResults)) {
    stop("plotResults must be a logical value")
  }

  seqLength <- length(seqCharacterVector)


  scoreVector <- KDNorm$V2[match(seqCharacterVector, KDNorm$V1)]

  hydropathyDF <- data.frame(Position = seq_len(seqLength),
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
#' Calculate the Mean Scaled Hydropathy
#'
#' This function utilizes the scaledHydropathyGlobal() function and
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
#' @examples
#' #Amino acid sequences can be character strings
#' aaString <- "ACDEFGHIKLMNPQRSTVWY"
#' #Amino acid sequences can also be character vectors
#' aaVector <- c("A", "C", "D", "E", "F",
#'               "G", "H", "I", "K", "L",
#'               "M", "N", "P", "Q", "R",
#'               "S", "T", "V", "W", "Y")
#'  #Alternativly, .fasta files can also be used by providing
#'
#' #Calculate the mean scaled hydropathy
#'  meanScaledHydropathy(aaString)
#'  meanScaledHydropathy(aaVector)

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
