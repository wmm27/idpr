sequenceMapCoordinates <-
  function(sequence,
           nbResidues = 30) {

    seqCharacterVector <- sequenceCheck(
      sequence = sequence,
      sequenceName = NA,
      method = "Stop",
      outputType = "Vector",
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
