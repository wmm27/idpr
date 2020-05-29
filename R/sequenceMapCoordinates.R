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
#'
#' @seealso \code{\link{sequenceMapCoordinates}} for mapping coordinates
#' @export
#'




sequenceMapCoordinates <-
  function(sequence,
           nbResidues = 30) {

    seqCharacterVector <- sequenceCheck(
      sequence = sequence,
      method = "stop",
      outputType = "vector",
      supressOutputMessage = T)

    seqLength <- length(seqCharacterVector)
    seqDF <- data.frame(Position = c(1:seqLength),
                        AA = seqCharacterVector)
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

    return(seqDF)
  }
