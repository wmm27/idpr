#' Protein Charge Calculation, Globally
#'
#' This function will determine the charge of a peptide using the
#'   Henderson-Hasselbalch Equation. The output is a data frame (default) or q
#'   plot of charge calculations along the peptide sequence. Charges are
#'   determined globally, or along the entire chain.
#'
#' @param sequence amino acid sequence as a character string or vector of
#'    individual residues. alternatively, a character string of the path to a
#'    .fasta / .fa file
#' @param pKaSet A character string or data frame. "IPC_protein" by default.
#'   Character string to load specific, preloaded pKa sets.
#'    c("EMBOSS", "DTASelect", "Solomons", "Sillero", "Rodwell",
#'     "Lehninger", "Toseland", "Thurlkill", "Nozaki", "Dawson",
#'     "Bjellqvist", "ProMoST", "Vollhardt", "IPC_protein", "IPC_peptide")
#'    Alternativly, the user may supply a custom pKa dataset.
#'    The format must be a data frame where:
#'    Column 1 must be a character vector of residues named "AA" AND
#'    Column 2 must be a numeric vector of pKa values.
#' @param pH numeric value, 7.0 by default.
#'   The environmental pH used to calculate residue charge.
#' @param plotResults logical value, FALSE by default.
#'   This determines what is returned. If \code{plotResults = FALSE}, a
#'   data frame is returned with the position, residue, and charge (-1 to +1).
#'   If  \code{plotResults = TRUE}, a graphical output is returned (ggplot)
#'   showing the charge distribution.
#' @param includeTermini,sumTermini Logical values, both TRUE by default. This
#'   determines how the calculation handles the N- and C- terminus.
#'   includeTermini determines if the calculation will use the charge of the
#'   amine and carboxyl groups at the ends of the peptide (When TRUE). These
#'   charges are ignored when \code{includeTermini = FALSE}. sumTermini
#'   determines if the charge of the first (likely Met, therefore uncharged),
#'   and final residue (varies) will be added to the termini charges, or if the
#'   N and C terminus will be returned as separate residues.
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
#' @return  If \code{plotResults = FALSE}, a data frame
#'   is returned with the position, residue, and charge (-1 to +1). If
#'   \code{plotResults = TRUE}, a graphical output is returned (ggplot) showing
#'   the charge distribution.
#' @family charge functions
#' @seealso \code{\link{pKaData}} for residue pKa values and
#'   \code{\link{hendersonHasselbalch}} for charge calculations.
#' @export
#' @examples
#'  #Amino acid sequences can be character strings
#' aaString <- "ACDEFGHIKLMNPQRSTVWY"
#' #Amino acid sequences can also be character vectors
#' aaVector <- c("A", "C", "D", "E", "F",
#'               "G", "H", "I", "K", "L",
#'               "M", "N", "P", "Q", "R",
#'               "S", "T", "V", "W", "Y")
#' #Alternativly, .fasta files can also be used by providing
#' #a character string of the path to the file.
#' exampleDF <- chargeCalculationGlobal(aaString)
#' head(exampleDF)
#' exampleDF <- chargeCalculationGlobal(aaVector)
#' head(exampleDF)
#'
#'
#' #Changing pKa set or pH used for calculations
#' exampleDF_pH5 <- chargeCalculationGlobal(aaString,
#'                                          pH = 5)
#' head(exampleDF_pH5)
#' exampleDF_pH7 <- chargeCalculationGlobal(aaString,
#'                                          pH = 7)
#' head(exampleDF_pH7)
#' exampleDF_EMBOSS <- chargeCalculationGlobal(aaString,
#'                                             pH = 7,
#'                                             pKa = "EMBOSS")
#' head(exampleDF_EMBOSS)
#'
#' #If the termini charge should not be included with includeTermini = F
#' exampleDF_NoTermini <- chargeCalculationGlobal(aaString,
#'                                                includeTermini = FALSE)
#' head(exampleDF_NoTermini)
#'
#' #and how the termini should be handeled with sumTermini
#' exampleDF_SumTermini <- chargeCalculationGlobal(aaString,
#'                                                 sumTermini = TRUE)
#' head(exampleDF_SumTermini)
#' exampleDF_SepTermini <- chargeCalculationGlobal(aaString,
#'                                                 sumTermini = FALSE)
#' head(exampleDF_SepTermini)
#'
#' #plotResults = TRUE will output a ggplot as a line plot
#'   chargeCalculationGlobal(aaString,
#'                           plot = TRUE)
#'
#'   #since it is a ggplot, you can change or annotate the plot
#'   gg <- chargeCalculationGlobal(aaVector,
#'                                 window = 3,
#'                                 plot = TRUE)
#'   gg <- gg + ggplot2::ylab("Residue Charge")
#'   gg <- gg + ggplot2::geom_text(data = exampleDF,
#'                                 ggplot2::aes(label = AA,
#'                                              y = Charge + 0.1))
#'   plot(gg)
#' #alternativly, you can pass the data frame to sequenceMap()
#' sequenceMap(sequence = exampleDF$AA,
#'             property = exampleDF$Charge)

chargeCalculationGlobal <- function(
  sequence,
  pKaSet = "IPC_protein",
  pH = 7.0,
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
    supressOutputMessage = TRUE)

  seqLength <- length(seqCharacterVector)
  positionVector <- seq_len(seqLength)

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
    if (dim(pKaSet)[2] >= 2) {
      pKaUsed <- pKaSet
      names(pKaUsed)[c(1, 2)] <- c("AA", "pKa")
    } else {
      stop("Custom pKaSet must be a data frame with 2 (or more) columns.
           Column 1 must be a character vector of residues.
           Column 2 must be a numeric vector of pKa values.")
    }
  }

  #------ Calculates the charge for each residue at the pH given
  iterations <- nrow(pKaUsed)
  pKaUsed$charge <- rep(NA, iterations)

  for (i in seq_len(iterations)) {

    residue.i <- pKaUsed[i, 1]
    pKa.i <- pKaUsed[i, 2]

    charge.i <- hendersonHasselbalch(
      pKa = as.numeric(pKa.i),
      pH = pH,
      residue = as.character(residue.i))

    pKaUsed$charge[i] <- charge.i
  }
  pKaUsed <- as.data.frame(pKaUsed)
  #----- matches the pKa set to the sequence
  pKaMatched <- pKaUsed[match(seqCharacterVector, pKaUsed[, 1]), 3]
  pKaMatched <- as.numeric(pKaMatched)
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
      plotTitle <- paste0("Distribution of Charged Residues in ", proteinName)
    } else {
      plotTitle <- "Distribution of Charged Residues"
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

#' Charge Calculation Along a Protein Sequence
#'
#' This calculates the charge, as determined by the Henderson-Hasselbalch
#'   equation, for each window along the sequence. This function uses a
#'   sliding window. The output is either a graph or a data frame of
#'   calculated charges.
#' @inheritParams chargeCalculationGlobal
#' @param sequence amino acid sequence as a single character string
#'   or vector of single characters.
#'   It also supports a single character string that specifies
#'   the location of a .fasta or .fa file.
#' @param window a positive, odd integer. 7 by default.
#'   Sets the size of sliding window, must be an odd number.
#'   The window determines the number of residues to be analyzed and averaged
#'   for each position along the sequence.
#' @param plotResults logical value. TRUE by default.
#'   If \code{plotResults = TRUE}, a ggplot of window charges are returned.
#'   If \code{plotResults = FALSE}, a data frame of window charges are returned.
#' @param proteinName character string, optional. Used to add protein name
#'   to the title in ggplot. Ignored if \code{plotResults = FALSE}.
#' @return see plotResults argument
#' @family charge functions
#' @seealso \code{\link{pKaData}} for residue pKa values and citations. See
#'   \code{\link{hendersonHasselbalch}} for charge calculations.
#' @export
#' @examples
#'  #Amino acid sequences can be character strings
#' aaString <- "ACDEFGHIKLMNPQRSTVWY"
#' #Amino acid sequences can also be character vectors
#' aaVector <- c("A", "C", "D", "E", "F",
#'               "G", "H", "I", "K", "L",
#'               "M", "N", "P", "Q", "R",
#'               "S", "T", "V", "W", "Y")
#' #Alternativly, .fasta files can also be used by providing
#' # a character string of the path to the file.
#' exampleDF <- chargeCalculationLocal(aaString)
#' exampleDF <- chargeCalculationLocal(aaVector)
#' head(exampleDF)
#'
#' #Changing window will alter the number of residues analyzed
#' exampleDF_window3 <- chargeCalculationLocal(aaString,
#'                                             window = 3)
#' head(exampleDF_window3)
#' exampleDF_window15 <- chargeCalculationLocal(aaString,
#'                                              window = 15)
#' head(exampleDF_window15)
#'
#' #Changing pKa set or pH used for calculations
#' exampleDF_pH5 <- chargeCalculationLocal(aaString,
#'                                         pH = 5)
#' head(exampleDF_pH5)
#' exampleDF_pH7 <- chargeCalculationLocal(aaString,
#'                                        pH = 7)
#' head(exampleDF_pH7)
#' exampleDF_EMBOSS <- chargeCalculationLocal(aaString,
#'                                            pH = 7,
#'                                            pKa = "EMBOSS")
#' head(exampleDF_EMBOSS)
#'
#' #plotResults = TRUE will output a ggplot
#'   chargeCalculationLocal(aaString,
#'                          plot = TRUE)
#'
#'   #since it is a ggplot, you can change or annotate the plot
#'   gg <- chargeCalculationLocal(aaVector,
#'                                window = 3,
#'                                plot = TRUE)
#'   gg <- gg + ggplot2::ylab("Local Charge")
#'   gg <- gg + ggplot2::geom_text(data = exampleDF_window3,
#'                                 ggplot2::aes(label = CenterResidue,
#'                                              y = windowCharge + 0.1))
#'  plot(gg)


chargeCalculationLocal <- function(
  sequence,
  window = 7,
  pH = 7.0,
  pKaSet = "IPC_protein",
  printCitation = FALSE,
  plotResults = FALSE,
  proteinName = NA,
  ...) {

  seqVector <- sequenceCheck(
    sequence = sequence,
    method = "stop",
    outputType = "vector",
    supressOutputMessage = TRUE)
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

  for (i in seq_len(iterations)) {
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
  for (i in seq_len(numberResiduesAnalyzed)) {
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


#' Protein Charge Calculation, Net Charge
#'
#' This function will determine the net charge of a peptide using the
#'   Henderson-Hasselbalch Equation. The output is a numeric value describing
#'   the total net charge or the average net charge.
#'
#' @inheritParams chargeCalculationGlobal
#' @param averaged logical value. FALSE by default.
#'   When \code{averaged = FALSE}, the total net charge is returned.
#'   When \code{averaged = TRUE}, the total net charge is averaged by the
#'   sequence length. This gives a value of -1 to +1.
#' @param includeTermini Logical value, TRUE by default. This
#'   determines how the calculation handles the N- and C- terminus.
#'   includeTermini determines if the calculation will use the charge of the
#'   amine and carboxyl groups at the ends of the peptide (When TRUE). These
#'   charges are ignored when \code{includeTermini = FALSE}.
#' @return numeric value. Either the net charge or average net charge, depending
#'   on the value of the averaged argument
#' @family charge functions
#' @seealso \code{\link{pKaData}} for residue pKa values and citations. See
#'   \code{\link{hendersonHasselbalch}} for charge calculations.
#' @export
#' @examples
#' #Amino acid sequences can be character strings
#' aaString <- "ACDEFGHIKLMNPQRSTVWY"
#' #Amino acid sequences can also be character vectors
#' aaVector <- c("A", "C", "D", "E", "F",
#'               "G", "H", "I", "K", "L",
#'               "M", "N", "P", "Q", "R",
#'               "S", "T", "V", "W", "Y")
#' #Alternativly, .fasta files can also be used by providing a character string
#'  # of the path to the file.
#'
#' #Calculate the Net Charge
#' netCharge(aaString,
#'           averaged = FALSE)
#' netCharge(aaVector,
#'           averaged = FALSE)
#'
#' #Calculate the Average Net Charge
#' netCharge(aaString,
#'           averaged = TRUE)
#' netCharge(aaVector,
#'           averaged = TRUE)
#'
#' #Change the pH
#' netCharge(aaString,
#'           pH = 8)
#' netCharge(aaString,
#'           pH = 7)
#' netCharge(aaString,
#'           pH = 5.5)
#'
#' #Specify which pKa set to use
#' netCharge(aaString,
#'           pKaSet = "IPC_protein") #Default
#' netCharge(aaString,
#'           pKaSet = "IPC_peptide")
#' netCharge(aaString,
#'           pKaSet = "Dawson")
#' netCharge(aaString,
#'           pKaSet = "EMBOSS")
#'
#' #Should the termini be included in charge calculations?
#' netCharge(aaString,
#'           includeTermini = TRUE) #Default
#' netCharge(aaString,
#'           includeTermini = FALSE)


netCharge <-
  function(sequence,
           pKaSet = "IPC_protein",
           pH = 7.0,
           includeTermini = TRUE,
           averaged = FALSE) {

    if (!is.logical(averaged)) {
     stop("averaged must be a logical value")
    }

    chargeDF <- chargeCalculationGlobal(
      sequence = sequence,
      pKaSet = pKaSet,
      pH = pH,
      plotOutput = FALSE,
      includeTermini = includeTermini,
      sumTermini = TRUE)

    netCharge <- sum(chargeDF$Charge)

    if (averaged) {
      seqLength <- nrow(chargeDF)
      avgNetCharge <- netCharge / seqLength
      return(avgNetCharge)
    } else {
      return(netCharge)
    }
  }
