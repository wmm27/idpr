#' Protein Charge Calculation, Globally
#'
#' This function will determine the charge of a peptide using the
#'   Henderson-Hasselbalch Equation. The output is a dataframe (default) or
#'   plot of charge calculations along the peptide sequence. Charges are
#'   determined globally, or along the entire chain.
#'
#' @param sequence amino acid sequence as a character string or vector of
#'    individual residues. alternativly a character string of the path to a
#'    .fasta / .fa file
#' @param pKaSet A character string or data frame. "IPC_protein" by default.
#'   Character string to load specific, preloaded pKa sets.
#'    c("AA", "EMBOSS", "DTASelect", "Solomons", "Sillero", "Rodwell",
#'     "Lehninger", "Toseland", "Thurlkill", "Nozaki", "Dawson",
#'     "Bjellqvist", "ProMoST", "IPC_protein", "IPC_peptide")
#'    Alternativly, the user may supply a custom pKa dataset.
#'    The format must be a data frame where:
#'    Column 1 must be a character vector of residues AND
#'    Column 2 must be a numeric vector of pKa values.
#' @param pH numeic value, 7.2 by default.
#'   The envioronmental pH used to calculate residue charge.
#' @param plotResults logical value, FALSE by default.
#'   This determines what is returned. If \code{plotResults = FALSE}, a
#'   dataframe is returned with the position, residue, and charge (-1 to +1). If
#'   \code{plotResults = TRUE}, a graphical output is returned (ggplot) showing
#'   the charge distribution.
#' @param includeTermini,sumTermini Logical values, both TRUE by default. This
#'   determines how the calculation handles the N- and C- terminus.
#'   includeTermini determines if the calculation will use the charge of the
#'   amine and carboxyl groups at the ends of the peptide (When TRUE). These
#'   charges are ignored when \code{includeTermini = FALSE}. sumTermini
#'   determines if the charge of the first (likely Met, therefore uncharged),
#'   and final residue (varies) will be added to the termini charges, or if the
#'   N and C terminus will be returned as seperate residues.
#'   When \code{sumTermini = TRUE}, charges are summed. When
#'   \code{sumTermini = FALSE}, the N and C terminus are added as a unique
#'   residue in the DF. This will impact averages by increasing the sequence
#'   length by 2. sumTermini is ignored if \code{includeTermini = FALSE}.
#' @param printCitation Logical value. FALSE by default.
#'   When \code{printCitation = TRUE} the citation for the pKa set is printed.
#'   This allows for the user to easily obtain the dataset citation.
#'   Will not print if there is a custom dataset.
#' @param proteinName character string with length = 1.
#'   optional setting to include the name in the plot title.
#' @param ... any additional parameters, especially those for plotting.
#' @return  If \code{plotResults = FALSE}, a dataframe
#'   is returned with the position, residue, and charge (-1 to +1). If
#'   \code{plotResults = TRUE}, a graphical output is returned (ggplot) showing
#'   the charge distribution.
#' @family charge functions
#' @seealso \code{\link{pKaData}} for residue pKa values and
#'   \code{\link{hendersonHasselbalch}} for charge calculations.
#' @export


chargeCalculationGlobal <- function(
sequence,
pKaSet = "IPC_protein",
pH = 7.2,
plotResults = FALSE,
includeTermini = TRUE,
sumTermini = TRUE,
proteinName = NA,
printCitation = FALSE,
...) {

  seqCharacterVector <- sequenceCheck(
    sequence = sequence,
    method = "stop",
    outputType = "vector",
    supressOutputMessage = T)

  seqLength <- length(seqCharacterVector)
  positionVector <- 1:seqLength

  if (!all(c(is.logical(plotResults),
             is.logical(includeTermini),
             is.logical(sumTermini),
             is.logical(printCitation)))) {
    stop(" plotResults, includeTermini, sumTermini, and printCitation
         must be logical values")
  }

  #----- to get pKa Set

  if (is.character(pKaSet) && length(pKaSet) == 1) {
    if (pKaSet %in% names(pKaData)) { #This is the section for supplied pKa sets
      setVector <- names(pKaData) %in% c("AA", pKaSet)
      pKaUsed <- pKaData[, setVector]

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
    if (dim(pKaUsed)[2] >= 2) {
      pKaUsed <- pKaSet
    } else {
      stop("Custom pKaSet must be a data frame with 2 (or more) columns.
           Column 1 must be a character vector of residues.
           Column 2 must be a numeric vector of pKa values.")
    }
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

  #----- matches the pKa set to the sequence
  pKaMatched <- pKaUsed[match(seqCharacterVector, pKaUsed$AA), 3]
  pKaMatched <- as.numeric(pKaMatched$charge)
  pKaMatched[is.na(pKaMatched)] <- 0

  if (includeTermini) {
    nTerminus <- c(pKaUsed[pKaUsed$AA == "NH2", 3],
                   pKaUsed[pKaUsed$AA == "NH3", 3])
    nTerminus <- as.numeric(nTerminus)
    nTerminus <- nTerminus[!is.na(nTerminus)][1] #incase two terminus values

    cTerminus <- c(pKaUsed[pKaUsed$AA == "COOH", 3],
                   pKaUsed[pKaUsed$AA == "COO", 3])
    cTerminus <- as.numeric(cTerminus)
    cTerminus <- cTerminus[!is.na(cTerminus)][1] #incase two terminus values

    if (sumTermini) {
      pKaMatched[1] <- pKaMatched[1] + nTerminus
      pKaMatched[seqLength] <- pKaMatched[seqLength] + cTerminus
    } else {
      pKaMatched <- c(nTerminus, pKaMatched, cTerminus)
      positionVector <- c(0 : (seqLength + 1))
      seqCharacterVector <- c("NH3", seqCharacterVector, "COO")
    }
  }

  chargeDF <- data.frame(Position = positionVector,
                         AA = seqCharacterVector,
                         Charge = pKaMatched)

  totalCharge <- sum(chargeDF$Charge)

  if (plotResults) {

    if (!is.na(proteinName)) {
      plotTitle <- paste0("Disribution of Charged Residues in ", proteinName)
    } else {
      plotTitle <- "Disribution of Charged Residues"
    }
    plotSubtitle <- paste0("Net Charge = ",
                          totalCharge)


    gg <-  sequencePlot(
      position = chargeDF$Position,
      property = chargeDF$Charge,
      hline = 0,
      dynamicColor = chargeDF$Charge,
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
