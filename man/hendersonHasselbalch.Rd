% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hendersonHasselbalch.R
\name{hendersonHasselbalch}
\alias{hendersonHasselbalch}
\title{Henderson-Hasselbalch Equation}
\usage{
hendersonHasselbalch(pKa, pH = 7, residue)
}
\arguments{
\item{pKa}{numeric value. The point where A- = HA.}

\item{pH}{numeric value. The pH of the environment. 7.0 by default}

\item{residue}{individual character or character string.
accepted values are the exact aa c("C", "D", "E", "H", "K", "R", "Y"),
termini c("COOH","COO","NH2","NH3"), or a
general property c("acid", "base", "negative", "positive").}
}
\value{
a numeric value giving the ratio of charged to uncharged residues.
}
\description{
This function calculates the ionic charge of a residue at a specific pH
  when given the pKa. \eqn{pH = pKa + log([A-] / [HA])}
  Known, charged residues are accepted as well as
  the protein termini and general property to allow customized calculations.
  The output is a ratio comparing acid to conjugate base for acidic residues
  or a ratio comparing conjugate base to acid for basic residues.
}
\examples{
#Calculating Lysine charge using the EMBOSS pKa data
EMBOSS_pKa <- pKaData[, 1:2]
EMBOSS_pKa

Lys_pKa <- EMBOSS_pKa[EMBOSS_pKa$AA == "K", ]
Lys_pKa$EMBOSS #This is Lysines pKa

hendersonHasselbalch(
  pKa = as.numeric(Lys_pKa$EMBOSS),
  pH = 7.0,
  residue = "K")

#residue = supports general properties as well
hendersonHasselbalch(
  pKa = as.numeric(Lys_pKa$EMBOSS),
  pH = 7.0,
  residue = "base")

hendersonHasselbalch(
  pKa = as.numeric(Lys_pKa$EMBOSS),
  pH = 7.0,
  residue = "positive")

#CALCULATIONS ARE DEPENDENT ON RESIDUE PROPERTY!
hendersonHasselbalch(
  pKa = as.numeric(Lys_pKa$EMBOSS),
  pH = 7.0,
  residue = "acid") #Inaccurate Description

#You can also calculate charge at different pHs
hendersonHasselbalch(
  pKa = as.numeric(Lys_pKa$EMBOSS),
  pH = 5.5,
  residue = "K")
hendersonHasselbalch(
  pKa = as.numeric(Lys_pKa$EMBOSS),
  pH = 8,
  residue = "K")

}
\seealso{
\code{\link{pKaData}} for residue pKa values and citations. See
 other charge functions for use.

Other charge functions: 
\code{\link{chargeCalculationGlobal}()},
\code{\link{chargeCalculationLocal}()},
\code{\link{netCharge}()}
}
\concept{charge functions}
