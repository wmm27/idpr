#' Henderson-Hasselbalch Equation
#'
#' This function calculates the ionic charge of a residue at a specific pH
#'   when given the pKa. \eqn{pH = pKa + log([A-] / [HA])}
#'   Known, charged residues are accepted as well as
#'   the protein termini and general property to allow customized calculations.
#'   The output is a ratio comparing acid to conjugate base for acidic residues
#'   or a ratio comparing conugate base to acid for basic residues.
#'
#' @param pKa numeric value. The point where A- = HA.
#' @param pH numeric value. The pH of the enviornment. 7.2 by default
#' @param residue individual character or character string.
#'   accepted values are the exact aa c("C", "D", "E", "H", "K", "R", "Y"),
#'   termini c("COOH","COO","NH2","NH3"), or a
#'   general property c("acid","base","negative", "positive").
#' @return a numeric value giving the ratio of charged to uncharged residues.
#' @family charge functions
#' @seealso \code{\link{pKaData}} for residue pKa values and citations. See
#'  other charge functions for use.
#' @export
#' @examples
#' #Calculating Lysine charge using the EMBOSS pKa data
#' EMBOSS_pKa <- pKaData[, 1:2]
#' EMBOSS_pKa
#'
#' Lys_pKa <- EMBOSS_pKa[EMBOSS_pKa$AA == "K", ]
#' Lys_pKa$EMBOSS #This is Lysines pKa
#'
#' hendersonHasselbalch(
#'   pKa = as.numeric(Lys_pKa$EMBOSS),
#'   pH = 7.2,
#'   residue = "K")
#'
#' #residue = supports general properties as well
#' hendersonHasselbalch(
#'   pKa = as.numeric(Lys_pKa$EMBOSS),
#'   pH = 7.2,
#'   residue = "base")
#'
#' hendersonHasselbalch(
#'   pKa = as.numeric(Lys_pKa$EMBOSS),
#'   pH = 7.2,
#'   residue = "positive")
#'
#' #CALCULATIONS ARE DEPENDENT ON RESIDUE PROPERTY!
#' hendersonHasselbalch(
#'   pKa = as.numeric(Lys_pKa$EMBOSS),
#'   pH = 7.2,
#'   residue = "acid") #Inaccurate Description
#'
#' #You can also calculate charge at different pHs
#' hendersonHasselbalch(
#'   pKa = as.numeric(Lys_pKa$EMBOSS),
#'   pH = 5.5,
#'   residue = "K")
#' hendersonHasselbalch(
#'   pKa = as.numeric(Lys_pKa$EMBOSS),
#'   pH = 8,
#'   residue = "K")
#'
hendersonHasselbalch <- function(
  pKa,
  pH = 7.2,
  residue) {
  #residue supports the exact aa c('C','D','E','H','K','R','Y'),
  #termini c("COOH","COO","NH2","NH3"), or
  #property c('acid','base','positive','negative')

  validAcidicResidues <- c("D", "E", "C", "Y",
                           "acid", "negative",
                           "COOH", "COO")
  validBasicResidues <-  c("H", "K", "R",
                           "base", "positive",
                           "NH2", "NH3")
  #---- Validating Input
  if (!residue %in% c(validAcidicResidues, validBasicResidues)) {
    stop("Please set residue equal to an accepted value.")
  }
  if (is.numeric(pKa) == FALSE ||
      is.numeric(pH) == FALSE) {
    stop("pKa and pH each require one numeric value")
  }

  #------ Calculations depending on residue type
  if (residue %in% validAcidicResidues) {
    charge <- -1 / (1 + 10 ^ (pKa - pH))
    return(charge)
  }
  if (residue %in% validBasicResidues) {
    charge <- 1 / (1 + 10 ^ (pH - pKa))
    return(charge)
  }
}
