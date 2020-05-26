#' Charge Calculation Along a Protein Sequence
#'
#' This calculates the charge, as determined by the Henderson-Hasselbalch
#'   equation, for each window along the sequence. This function uses a
#'   sliding window. The output is either a graph or dataframe of calculations.

#' @inheritParams chargeCalculationGlobal
#'
#' @param sequence amino acid sequence as a single character string
#'   or vector of single characters.
#'   It also supports a single character string that specifies
#'   the locaion of a .fasta or .fa file.
#' @param window a positive, odd integer. 7 by default.
#'   Sets the size of sliding window, must be an odd number.
#'   The window determines the number of residues to be analyzed and averaged
#'   for each position along the sequence.
#' @param plotResults logical value. TRUE by default.
#'   If \code{plotResults = TRUE}, a ggplot of local charge is returned
#'   If \code{plotResults = FALSE}, a dataframe of window charges are returned.
#' @param proteinName character string, optional. Used to add protein name
#'   to the title in ggplot. Ignored if \code{plotResults = FALSE}.
#' @return see plotResults argument
#' @family charge functions
#' @seealso \code{\link{pKaData}} for residue pKa values and citations. See
#'   \code{\link{hendersonHasselbalch}} for charge calculations.
#' @export



chargeCalculationLocal <- function(
  sequence,
  window = 7,
  pH = 7.2,
  pKaSet = "IPC_protein",
  printCitation = FALSE,
  plotResults = FALSE,
  proteinName = NA,
  ...) {

  seqVector <- sequenceCheck(
    sequence = sequence,
    sequenceName = proteinName,
    method = "Stop",
    outputType = "Vector",
    supressOutputMessage = T)
  seqLength <- length(seqVector)


  if ((window %% 2) == 0) {
    stop("Window must be an odd number")
  }

  if (!all(c(is.logical(plotResults),
      is.logical(printCitation)))) {
    stop("plotResults and printCitation must be logical values")
  }

  if (is.character(pKaSet) && length(pKaSet) == 1) {
    if (pKaSet %in% names(pKaData)) { #This is the section for supplied pKa sets
      setVector <- names(pKaData) %in% c("AA", pKaSet)
      pKaUsed <- pKaData[, setVector]
      pKaUsed <- as.data.frame(pKaUsed)
      if (printCitation) {
        pKaCitation <- pKaUsed[pKaUsed$AA == "citation", 2]
        print(pKaCitation)
      }

      pKaUsed <- pKaUsed[!pKaUsed$AA == "citation", ]

    } else {
      stop("Invalid pKa Set Provided")
    }
  }

  if (is.data.frame(pKaSet)) {
      pKaUsed <- pKaSet
   }


  #------ Calculates the charge for each residue at the pH given
  iterations <- nrow(pKaUsed)
  pKaUsed$charge <- rep(NA, iterations)

  for (i in 1:iterations) {
    residue.i <- pKaUsed[i, 1]
    pKa.i <- pKaUsed[i, 2]

    charge.i <- hendersonHasselbalch(
      pKa = as.numeric(pKa.i),
      pH = pH,
      residue = as.character(residue.i))

    pKaUsed$charge[i] <- charge.i
  }
   pKaUsed <- as.data.frame(pKaUsed)
  #the number of possible windows
  numberResiduesAnalyzed <- seqLength - (window - 1)

  #empty vectors to loop into
  windowVector <- rep(NA, numberResiduesAnalyzed)
  scoreVector <- rep(NA, numberResiduesAnalyzed)

  #----- Analysis
  for (i in 1:numberResiduesAnalyzed) {
    windowBegining <- i
    windowEnd      <- i + window - 1
    windowResidues <- seqVector[windowBegining:windowEnd]
    windowVector[i] <- paste(windowResidues, collapse = "")

    pKaMatched <- pKaUsed[match(windowResidues, (pKaUsed[, 1])), 3]
    pKaMatched <- as.numeric(pKaMatched)
    pKaMatched[is.na(pKaMatched)] <- 0

    scoreVector[i] <- sum(pKaMatched) / window
  }
  numberVector <- ((window + 1) / 2) : (seqLength - (window - 1) / 2)
  residueVector <- seqVector[numberVector]

  chargeDF <- data.frame(Position = numberVector,
                         CenterResidue = residueVector,
                         Window = windowVector,
                         windowCharge = scoreVector)

  if (plotResults) {
    if (!is.na(proteinName)) {
      plotTitle <- paste0("Calculation of Local Charge in ", proteinName)
    } else {
      plotTitle <- "Calculation of Local Charge"
    }
    netChargeValue <- netCharge(sequence = seqVector)
    netChargeValue <- round(netChargeValue, 3)
    plotSubtitle <- paste0("Window Size = ",
                           window,
                           " ; Net Charge = ",
                           netChargeValue)
    gg <-  sequencePlot(
      position = chargeDF$Position,
      property = chargeDF$windowCharge,
      hline = 0,
      dynamicColor = chargeDF$windowCharge,
      customColors = c("blue", "red", "grey65"),
      customTitle = NA,
      propertyLimits = c(-1, 1))

    gg <- gg + ggplot2::labs(title = plotTitle,
                             subtitle = plotSubtitle)
    return(gg)
  } else {
    return(chargeDF)
  }
}
