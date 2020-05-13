#' Protein Charge Calculation, Net Charge
#'
#' This function will determine the net charge of a peptide using the
#'   Henderson-Hasselbalch Equation. The output is a numeric value. The output
#'   is the total net charge, but can be the average net charge.
#'
#' @inheritParams chargeCalculationGlobal
#' @param averaged logical value. FALSE by default.
#'   When \code{averaged = FALSE}, the total net charge is returned.
#'   When \code{averaged = TRUE}, the total net charge is averaged by the
#'   sequence length. This gives a value of -1 to +1.
#' @return numeric value. Either the net charge or average net charge, depending
#'   on the value of the averaged argument
#' @family charge functions
#' @seealso \code{\link{pKaData}} for residue pKa values and citations. See
#'   \code{\link{hendersonHasselbalch}} for charge calculations.
#' @export


netCharge <-
  function(sequence,
           pKaSet = "IPC_protein",
           pH = 7.2,
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
