#' Sequence Map Coordinates
#'
#' This is a function used to create a coordinate grid for the
#'   \code{\link{sequenceMap}} function. It is based on the length of the
#'   sequence being mapped, and how many residues per line are specified.
#'   The function wraps the sequence to have a number of columns that is
#'   the sequence length / number of residues per row, rounded up. \cr\cr
#'   This is intended for use within the sequenceMap function, however, this
#'   can also be used to identify the coordinates of residues within the ggplot
#'   coordinate plane for addition annotations.
#'
#' @inheritParams sequenceMap
#' @inheritParams sequenceCheck
#'
#' @return A data frame with rows containing the amino acid sequence, residue
#'   position within the sequence, as well as the row and column of each
#'   residue within the ggplot output of sequenceMap().
#' @seealso \code{\link{sequenceMapCoordinates}} for mapping coordinates
#' @export
#' @examples
#' #Amino acid sequences can be character strings
#' aaString <- "ACDEFGHIKLMNPQRSTVWY"
#' #Amino acid sequences can also be character vectors
#' aaVector <- c("A", "C", "D", "E", "F",
#'            "G", "H", "I", "K", "L",
#'            "M", "N", "P", "Q", "R",
#'            "S", "T", "V", "W", "Y")
#' #Alternativly, .fasta files can also be used by providing
#' ##The path to the file as a character string
#'
#' exampleDF <- sequenceMapCoordinates(aaString,
#'                                   nbResidues = 10)
#' head(exampleDF)
#'
#' exampleDF <- sequenceMapCoordinates(aaVector,
#'                                  nbResidues = 10)
#' head(exampleDF)
#'
#' #Getting a data frame for plotting with sequenceMap
#'  tendencyDF <- structuralTendency(sequence = aaVector)
#'  #Making a sequenceMap ggplot to annotate
#' gg <- sequenceMap(sequence = tendencyDF$AA,
#'             property = tendencyDF$Tendency,
#'              nbResidues = 3,
#'               labelType = "both")
#'
#' #Change the nbResidues to correspond to the sequenceMap setting
#' mapCoordDF <- sequenceMapCoordinates(aaVector,
#'                                   nbResidues = 3)
#' head(mapCoordDF)
#'
#' #subsetting for positive residues
#' mapCoordDF_subset <- mapCoordDF$AA %in% c("K", "R", "H")
#' mapCoordDF_subset <- mapCoordDF[mapCoordDF_subset,]
#'
#' #use mapCoordDF to annotate positive residues with a plus
#' library(ggplot2)
#' gg <- gg + geom_point(inherit.aes = FALSE,
#'                     data = mapCoordDF_subset,
#'                    aes(x = col + 0.5, #to center on the residue
#'                        y = row + 0.2), #to move above on the residue
#'                    color = "purple",
#'                    size = 3,
#'                    shape = 3)
#' plot(gg)

sequenceMapCoordinates <-
  function(sequence,
           nbResidues = 30) {

    seqCharacterVector <- sequenceCheck(
      sequence = sequence,
      method = "stop",
      outputType = "vector",
      supressOutputMessage = TRUE)

    seqLength <- length(seqCharacterVector)
    seqDF <- data.frame(Position = seq_len(seqLength),
                        AA = seqCharacterVector)
    nRows <- ceiling(seqLength / nbResidues)
    rowVector <- rep(1, seqLength)
    colVector <- seq_len(seqLength)
    if (nRows > 1) {
      for (i in seq_len(nRows)) {
        iteration <- i - 1
        row.i <- nRows - (iteration)
        start.i <- 1 + nbResidues * (iteration)
        end.i <- nbResidues +  nbResidues * (iteration)
        rowVector[start.i:end.i] <- row.i
        colVector[start.i:end.i] <- seq_len(nbResidues)
      }
    }

    seqDF$row <- rowVector[seq_len(seqLength)]
    seqDF$col <- colVector[seq_len(seqLength)]

    return(seqDF)
  }


#' Sequence Map Function
#'
#' This is a graphical function used to visualize data along an
#'   amino acid sequence.
#' @importFrom ggplot2 ggplot aes geom_bin2d theme geom_text aes_
#'
#' @inheritParams sequenceCheck
#' @param property a vector with length equal to sequence length.
#'   This is what is visualized on the function.
#'   Can be discrete or continuous values.
#' @param nbResidues numeric value, 30 by default.
#'   The number of residues to display on each row of the plot.
#'   It is not recommended to be over 50 or under 10 for standard sequences.
#'   Optimal value may vary between sequences of extreme lengths.
#' @param labelType character string, "both" by default.
#'   accepted values are \code{labelType = c("both", "AA", "number", "none")}.
#'   "both" shows both amino acid residue and residue number. "AA" and "number"
#'   show either the amino acid residue or the residue number, respectively.
#'   "none" only shows graphical values without labels.
#'   NOTE: When using "both", *everyN*, *labelLocation*, and *rotationAngle*
#'   all require vectors of length = 2 where the first value applies to the
#'   "AA" parameter and the second value applies to the "number" parameter.
#'   When using "AA" or "number, *everyN*, *labelLocation*, and *rotationAngle*
#'   require a single value. If a vector us provided, only the first value
#'   will be used.
#' @param everyN numeric value or vector of numeric values with length = 2.
#'   This is used to show every Nth amino acid and/or residue number.
#'   To show every value, set \code{everyN = 1} or \code{everyN = c(1, 1)}.
#' @param labelLocation character string or vector of character strings
#'   with length = 2. When \code{labelLocation = "on"}, the text is layered
#'   on top of the graphical output. When \code{labelLocation = "below"}, the
#'   text is placed below the graphical output.
#'   If \code{labelType = "both"}, do not set
#'   \code{labelLocation = c("on", "on")} or
#'   \code{labelLocation = c("below", "below")}.
#' @param rotationAngle numeric value or vector of numeric values
#'   with length = 2. This value is used to rotate text. Especially
#'   useful when printing many residue numbers.
#' @param customColors vector of colors as character strings. NA by default.
#' Used to support custom plot colors. If property is a discrete scale, a
#'   character vector of colors with length = number of unique discrete
#'   observations is required. If property is a continuous scale, a character
#'   vector of the colors for c("highColor","lowColor","midColor").
#'   Set NA to skip custom colors.
#'
#' @return A ggplot.
#' @export
#' @seealso \code{\link{sequenceMapCoordinates}} for mapping coordinates
#' @examples
#' #Get a data frame returned from another function
#' aaVector <- c("A", "C", "D", "E", "F",
#'               "G", "H", "I", "K", "L",
#'               "M", "N", "P", "Q", "R",
#'               "S", "T", "V", "W", "Y")
#' ## As a continuous property
#' exampleDF_cont <- chargeCalculationGlobal(sequence = aaVector)
#' head(exampleDF_cont)
#' ## Or as a discrete property
#' exampleDF_disc <- structuralTendency(sequence = aaVector)
#' head(exampleDF_disc)
#' sequenceMap(sequence = exampleDF_cont$AA,
#'           property = exampleDF_cont$Charge,
#'          nbResidues = 3,
#'          labelType = "both")
#'
#' sequenceMap(sequence = exampleDF_disc$AA,
#'          property = exampleDF_disc$Tendency,
#'          nbResidues = 3,
#'          labelType = "both")
#'
#' #Change the layout of labels
#' sequenceMap(
#' sequence = exampleDF_disc$AA,
#' property = exampleDF_disc$Tendency,
#' nbResidues = 3,
#' labelType = "AA") #Only AA residue Labels
#'
#' sequenceMap(
#' sequence = exampleDF_disc$AA,
#' property = exampleDF_disc$Tendency,
#' nbResidues = 3,
#' labelType = "number") #Only residue numner labels
#'
#' sequenceMap(
#' sequence = exampleDF_disc$AA,
#' property = exampleDF_disc$Tendency,
#' nbResidues = 3,
#' labelType = "none") #No labels
#'
#' #The text can also be rotated for ease of reading,
#'  ## espeically helpful for larger sequences.
#' sequenceMap(
#' sequence = exampleDF_disc$AA,
#' property = exampleDF_disc$Tendency,
#' labelType = "number",
#' labelLocation = "on",
#'   rotationAngle = 90)
#'
#' #Specify colors for continuous values
#'
#' sequenceMap(
#' sequence = exampleDF_cont$AA,
#' property = exampleDF_cont$Charge,
#' customColors = c("purple", "pink", "grey90"))
#'
#' #or discrete values
#' sequenceMap(
#' sequence = exampleDF_disc$AA,
#' property = exampleDF_disc$Tendency,
#' customColors = c("#999999", "#E69F00", "#56B4E9"))
#'
#'
#' #change the number of residues on each line with nbResidue
#' #or discrete values
#' sequenceMap(
#' sequence = exampleDF_disc$AA,
#' property = exampleDF_disc$Tendency,
#' nbResidues = 1)
#' sequenceMap(
#' sequence = exampleDF_disc$AA,
#' property = exampleDF_disc$Tendency,
#' nbResidues = 3)
#' sequenceMap(
#' sequence = exampleDF_disc$AA,
#' property = exampleDF_disc$Tendency,
#' nbResidues = 10)
#'
#'
#' #Use sequenceMapCoordinates for additional annotations
#' gg <- sequenceMap(sequence = exampleDF_disc$AA,
#'                property = exampleDF_disc$Tendency,
#'                nbResidues = 3,
#'                labelType = "both")
#'
#' #Change the nbResidues to correspond to the sequenceMap setting
#' mapCoordDF <- sequenceMapCoordinates(aaVector,
#'                                   nbResidues = 3)
#' head(mapCoordDF)
#'
#' #subsetting for positive residues
#' mapCoordDF_subset <- mapCoordDF$AA %in% c("K", "R", "H")
#' mapCoordDF_subset <- mapCoordDF[mapCoordDF_subset,]
#'
#use mapCoordDF to annotate positive residues with a plus
#' library(ggplot2)
#' gg <- gg + geom_point(inherit.aes = FALSE,
#'                     data = mapCoordDF_subset,
#'                    aes(x = col + 0.5, #to center on the residue
#'                        y = row + 0.2), #to move above on the residue
#'                    color = "purple",
#'                    size = 3,
#'                    shape = 3)
#' plot(gg)
#'

sequenceMap <- function(
  sequence,
  property,
  nbResidues = 30,
  labelType = "both",
  everyN = c(1, 10),
  labelLocation = c("on", "below"),
  rotationAngle = c(0, 0),
  customColors = NA) {

  seqDF <- sequenceMapCoordinates(sequence = sequence,
                                  nbResidues = nbResidues)
  seqLength <- nrow(seqDF)
  seqDF$Property <- property
  nRows <- ceiling(seqLength / nbResidues)

  # ---- plot
  gg <- ggplot(data = seqDF,
               aes_(x = ~ col,
                   y = ~ row,
                   fill = ~ Property)) +
    geom_bin2d(binwidth = c(0.99, 0.5),
               aes_(group = ~ Property)) +
    ggplot2::theme_void()  +
    ggplot2::ylim(0, nRows + 0.25) +
    theme(legend.position = "top")
  #--- labeling
  if (!labelType == "none") {
    if (labelType == "both") {
      #-- Adding AA for "both"
      labelingVector <- seqDF$AA
      #---- Section to alter labels to be every Nth residue
      includeVector <- c(TRUE, rep(FALSE, everyN[1] - 1))
      includeVectorLength <- length(includeVector)
      includeVectorRep <- seqLength / includeVectorLength
      includeVector <- rep(includeVector, includeVectorRep)
      labelingVector[includeVector == FALSE] <- NA
      #---
      if (labelLocation[1] == "on") {
        yAdjust <- -0.25
      }
      if (labelLocation[1] == "below") {
        yAdjust <- -0.70
      }
      gg <- gg + geom_text(aes(x = (col + 0.5) * 0.99,
                               y = row,
                               label = as.character(labelingVector),
                               angle = rotationAngle[1]),
                           size = 3,
                           nudge_y = yAdjust,
                           na.rm = TRUE)
      #-- Adding numbers for "both"
      labelingVectorNumbers <- seqDF$Position
      #---- Section to alter labels to be every Nth residue
      includeVectorNumber <- c(TRUE, rep(FALSE, everyN[2] - 1))
      includeVectorNumberLength <- length(includeVectorNumber)
      includeVectorRep <- seqLength / includeVectorNumberLength
      includeVectorNumber <- rep(includeVectorNumber, includeVectorRep)
      labelingVectorNumbers[includeVectorNumber == FALSE] <- NA
      #-----
      if (labelLocation[2] == "on") {
        yAdjust <- -0.25
      }
      if (labelLocation[2] == "below") {
        yAdjust <- -0.70
      }
      gg <- gg + geom_text(aes(x = (col + 0.5) * 0.99,
                               y = row,
                               label = as.character(labelingVectorNumbers),
                               angle = rotationAngle[2]),
                           size = 3,
                           nudge_y = yAdjust,
                           na.rm = TRUE)
    } else {
      if (labelType == "AA") {
        labelingVector <- seqDF$AA
      }
      if (labelType == "number") {
        labelingVector <- seqDF$Position
      }
      if (labelLocation[1] == "on") {
        yAdjust <- -0.25
      }
      if (labelLocation[1] == "below") {
        yAdjust <- -0.70
      }
      #---- Section to alter labels to be every Nth residue
      includeVector <- c(TRUE, rep(FALSE, everyN[1] - 1))
      includeVectorLength <- length(includeVector)
      includeVectorRep <- seqLength / includeVectorLength
      includeVector <- rep(includeVector, includeVectorRep)
      labelingVector[includeVector == FALSE] <- NA
      gg <- gg + geom_text(aes(x = (col + 0.5) * 0.99,
                               y = row,
                               label = as.character(labelingVector),
                               angle = rotationAngle[1]),
                           size = 3,
                           nudge_y = yAdjust,
                           na.rm = TRUE)
    }
  }
  if (!is.na(customColors[1])) {
    if (plyr::is.discrete(seqDF$Property)) {
      gg <- gg + ggplot2::scale_fill_manual(values = customColors)
    } else {
      gg <- gg + ggplot2::scale_fill_gradient2(high = customColors[1],
                                               low = customColors[2],
                                               mid = customColors[3])
    }
  }
  return(gg)
}
