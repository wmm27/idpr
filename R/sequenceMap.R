#' Sequence Map Function
#'
#' This is a graphical function used to visualize data along an
#'   amino acid sequence.
#' @importFrom ggplot2 ggplot aes geom_bin2d theme geom_text
#'
#' @inheritParams sequenceCheck
#' @param property a vector with length equal to sequence length.
#'   This is what is visualized on the function.
#'   Can be discrete or contnious values.
#' @param nbResidues numeric value, 30 by default.
#'   The number of residues to display on each row of the plot.
#'   It is not reccomended to be over 50 or under 10 for standard sequences.
#'   Optimal value may vary between sequences of extreme lengths.
#' @param labelType character string, "both" by default.
#'   accepted values are \code{labelType = c("both", "AA", "number", "none")}.
#'   "both" shows both amino acid residue and residue number. "AA" and "number"
#'   show either the amino acid residue or the residue number, respectivly.
#'   "none" only shows graphical values without labels.
#'   NOTE: When using "both", *everyN*, *labelLocation*, and *rotationAngle*
#'   all require vectors of length = 2 where the first value applies to the
#'   "AA" parameter and the second value applies to the "number" parameter.
#'   When using "AA" or "number, *everyN*, *labelLocation*, and *rotationAngle*
#'   require a single value. If a vector us provided, only the first value
#'   will be used.
#' @param everyN numeric value or vector of numeric values with length = 2.
#'   This is used to show evey Nth amino acid and/or residue number.
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
#'   obervations is required. If property is a continuous scale, a character
#'   vector of the colors for c("highColor","lowColor","midColor").
#'   Set NA to skip custom colors.
#'
#' @return A ggplot. Can be plotted immediatly or
#'   additional layers can be added.
#' @export
#'


sequenceMap <- function(
  sequence,
  property,
  nbResidues = 30, # number of residues per line on the plot
  labelType = "both",
  everyN = c(1, 10),
  labelLocation = c("on", "below"),
  rotationAngle = c(0, 0),
  customColors = NA) {

  seqCharacterVector <- sequenceCheck(
    sequence = sequence,
    sequenceName = NA,
    method = "Stop",
    outputType = "Vector",
    supressOutputMessage = T)
  seqLength <- length(seqCharacterVector)
  seqDF <- data.frame(Position = c(1:seqLength),
                      AA = seqCharacterVector,
                      Property = property)
  nRows <- ceiling(seqLength / nbResidues)
  rowVector <- rep(1, seqLength)
  colVector <- c(1:seqLength)
  if (nRows > 1) {
    for (i in 1:(nRows)) {
      iteration <- i - 1
      row.i <- nRows - (iteration)
      start.i <- 1 + nbResidues * (iteration)
      end.i <- nbResidues +  nbResidues * (iteration)
      rowVector[start.i:end.i] <- row.i
      colVector[start.i:end.i] <- c(1:nbResidues)
    }
  }

  seqDF$row <- rowVector[1:seqLength]
  seqDF$col <- colVector[1:seqLength]

  # ---- plot
  gg <- ggplot(data = seqDF,
               aes(x = col,
                   y = row,
                   fill = Property)) +
    geom_bin2d(binwidth = c(0.99, 0.5),
               aes(group = Property)) +
    ggplot2::theme_void()  +
   ggplot2::ylim(0, nRows + 0.25) +
    theme(legend.position = "top")
  #--- labeling
  if (!labelType == "none") {
    if (labelType == "both") {
      #-- Adding AA for "both"
      labelingVector <- seqDF$AA
      #---- Section to alter labels to be every Nth residue
      includeVector <- c(T, rep(F, everyN[1] - 1))
      includeVectorLength <- length(includeVector)
      includeVectorRep <- seqLength / includeVectorLength
      includeVector <- rep(includeVector, includeVectorRep)
      labelingVector[includeVector == F] <- NA
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
      includeVectorNumber <- c(T, rep(F, everyN[2] - 1))
      includeVectorNumberLength <- length(includeVectorNumber)
      includeVectorRep <- seqLength / includeVectorNumberLength
      includeVectorNumber <- rep(includeVectorNumber, includeVectorRep)
      labelingVectorNumbers[includeVectorNumber == F] <- NA
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
      includeVector <- c(T, rep(F, everyN[1] - 1))
      includeVectorLength <- length(includeVector)
      includeVectorRep <- seqLength / includeVectorLength
      includeVector <- rep(includeVector, includeVectorRep)
      labelingVector[includeVector == F] <- NA
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
