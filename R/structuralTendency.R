#' Structural Tendency of Amino Acid Residues
#'
#' Each amino acid residue has a  tendency to impact the order / disorder
#'   of the amino acid sequence. Some residues are disorder promoting, meaning
#'   they tend to favor disorder over ordered structures. These are typically
#'   hydophillic, charged, or small residues. Order promoting residues tend
#'   to be aliphatic, hydrophobic, aromatic, or form tertiary structures.
#'   Disorder neutral residues neither favor order or disordered strucutres.
#'
#' @inheritParams sequenceCheck
#' @param printCitation logical, FALSE by default.
#'    When \code{printCitation = TRUE}, a citation to the paper
#'    categorizing the strucutral impact of each residue is printed.
#' @param disorderPromoting,disorderNeutral,orderPromoting character vectors
#'    of individual residues to be matched with the input sequence. Defaults:
#'    \itemize{
#'      \item disorderPromoting = c("P", "E", "S", "Q", "K", "A", "G")
#'      \item orderPromoting =
#'         c("M", "N", "V", "H", "L", "F", "Y", "I", "W", "C")
#'      \item disorderNeutral = c("D", "T", "R")
#'    }
#'    It is not reccomended to change these.
#' @return a dataframe containing each residue from the sequence
#'   matched with its structural tendancy, defined by disorderPromoting,
#'   disorderNeutral, and orderPromoting.
#'   For convenient plotting see \code{\link{structuralTendencyPlot}}.
#' @family structural tendency
#' @references Kulkarni, Prakash, and Vladimir N. Uversky. "Intrinsically
#'   disordered proteins: the dark horse of the dark proteome."
#'   Proteomics 18.21-22 (2018): 1800061.
#'   \url{https://doi.org/10.1002/pmic.201800061}.
#' @export

structuralTendency <- function(
  sequence,
  disorderPromoting = c("P", "E", "S", "Q", "K", "A", "G"),
  disorderNeutral = c("D", "T", "R"),
  orderPromoting = c("M", "N", "V", "H", "L", "F", "Y", "I", "W", "C"),
  printCitation = F) {
  #-----
  seqCharacterVector <- sequenceCheck(
    sequence = sequence,
    method = "stop",
    outputType = "vector",
    supressOutputMessage = T)
  sequenceLength <- length(seqCharacterVector)

  #----- Matches residue with tendency
  structuralTendencyVector <- rep(NA, sequenceLength)

  disorderedResidues <- seqCharacterVector %in% disorderPromoting
  structuralTendencyVector[disorderedResidues] <- "Disorder Promoting"

  orderedResidues <- seqCharacterVector %in% orderPromoting
  structuralTendencyVector[orderedResidues] <- "Order Promoting"

  neutralResidues <- seqCharacterVector %in% disorderNeutral
  structuralTendencyVector[neutralResidues] <- "Disorder Neutral"

  #----- makes the data frame for output
  structureTendencyDF <- data.frame(Position = 1:sequenceLength,
                                    AA = seqCharacterVector,
                                    Tendency = structuralTendencyVector)

  structureTendencyDF$AA <- as.character(structureTendencyDF$AA)
  structureTendencyDF$Tendency <- as.character(structureTendencyDF$Tendency)
  structureTendencyDF$Position <- as.numeric(structureTendencyDF$Position)

  if (printCitation) {
    residueCitation <- "Kulkarni, P., & Uversky, V. N. (2018).
        Intrinsically disordered proteins: the dark horse of the dark proteome.
        Proteomics, 18(21-22), 1800061."
    print(residueCitation)
  }
  return(structureTendencyDF)
}
