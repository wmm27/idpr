#' Structural Tendency of Amino Acid Residues
#'
#' Each amino acid residue has a tendency to impact the order / disorder
#'   of the amino acid sequence. Some residues are disorder promoting, meaning
#'   they tend to favor disorder over ordered structures. These are typically
#'   hydrophilic, charged, or small residues. Order promoting residues tend
#'   to be aliphatic, hydrophobic, aromatic, or form tertiary structures.
#'   Disorder neutral residues neither favor order or disordered structures.
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
#'    It is not recommended to change these. These definitions are from
#'    Uversky (2013).
#' @return a data frame containing each residue from the sequence
#'   matched with its structural tendancy, defined by disorderPromoting,
#'   disorderNeutral, and orderPromoting.
#'   For convenient plotting see \code{\link{structuralTendencyPlot}}.
#' @family structural tendency
#' @references
#'   Uversky, V. N. (2013). A decade and a half of protein intrinsic disorder:
#'   Biology still waits for physics. Protein Science, 22(6), 693-724.
#'   \url{https://doi.org/10.1002/pro.2261}. \cr
#'   Kulkarni, Prakash, and Vladimir N. Uversky. "Intrinsically
#'   disordered proteins: the dark horse of the dark proteome."
#'   Proteomics 18.21-22 (2018): 1800061.
#'   \url{https://doi.org/10.1002/pmic.201800061}.
#' @export
#' @examples
#' #Amino acid sequences can be character strings
#' aaString <- "ACDEFGHIKLMNPQRSTVWY"
#' #Amino acid sequences can also be character vectors
#' aaVector <- c("A", "C", "D", "E", "F",
#'               "G", "H", "I", "K", "L",
#'               "M", "N", "P", "Q", "R",
#'              "S", "T", "V", "W", "Y")
#' #Alternativly, .fasta files can also be used by providing
#' ##The path to the file as a character string
#'
#' exampleDF <- structuralTendency(aaString)
#' head(exampleDF)
#' exampleDF <- structuralTendency(aaVector)
#' head(exampleDF)
#'
#' #If using a different definition of disordered residues
#' ##These residues are labeled as such from Dunker et al (2001),
#' ##"Intrinsically disordered protein."
#' exampleDF <- structuralTendency(aaString,
#'                disorderPromoting = c("A", "R", "G", "Q", "S", "P", "E", "K"),
#'                disorderNeutral = c("H", "M", "T", "D"),
#'                orderPromoting = c("W", "C", "F", "I", "Y", "V", "L", "N"))
#' head(exampleDF)

structuralTendency <- function(
  sequence,
  disorderPromoting = c("P", "E", "S", "Q", "K", "A", "G"),
  disorderNeutral = c("D", "T", "R"),
  orderPromoting = c("M", "N", "V", "H", "L", "F", "Y", "I", "W", "C"),
  printCitation = FALSE) {
  #-----
  seqCharacterVector <- sequenceCheck(
    sequence = sequence,
    method = "stop",
    outputType = "vector",
    supressOutputMessage = TRUE)
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
  structureTendencyDF <- data.frame(Position = seq_len(sequenceLength),
                                    AA = seqCharacterVector,
                                    Tendency = structuralTendencyVector)

  structureTendencyDF$AA <- as.character(structureTendencyDF$AA)
  structureTendencyDF$Tendency <- as.character(structureTendencyDF$Tendency)
  structureTendencyDF$Position <- as.numeric(structureTendencyDF$Position)

  if (printCitation) {
    residueCitation <- "Uversky, V. N. (2013).
     A decade and a half of protein intrinsic disorder:
     Biology still waits for physics.
     Protein Science, 22(6), 693-724.
     doi:10.1002/pro.2261"
    print(residueCitation)
  }
  return(structureTendencyDF)
}



#' Plotting Structural Tendency of Amino Acid Sequence
#'
#' Convenient graphing for the \code{\link{structuralTendency}} function.
#'
#' @param sequence amino acid sequence (or pathway to a fasta file)
#'    as a character string. Supports multiple sequences / files, as a
#'    character vector of strings.
#' @param graphType character string, required.
#'   graphType must be set to c("pie", "bar", "none").
#'   When \code{graphType = "pie"}, the output is a pie chart.
#'   When \code{graphType = "bar"}, the output is a bar chart.
#'   When \code{graphType = "none"}, the output is the data frame that would
#'   otherwise be used to plot the data.
#' @param summarize logical value, FALSE by default.
#'   When \code{summarize = TRUE}, each residue is aggregated into Disorder
#'   Tendency Groups. (See \code{\link{structuralTendency}} for more details).
#'   When \code{summarize = FALSE}, residue identity is preserved and
#'   the output is colored by Disorder Tendency Groups.
#' @param alphabetical logical value, FALSE by default.
#'   Order of residues on plot axis. Only relevant when
#'    \code{summarize = FALSE}, otherwise is ignored.
#'    If FALSE, ordering is grouped by Disorder Tendency (P, E, S, ..., W, C).
#'    If TRUE, the residues are ordered alphabetically (A, C, D, E, ..., W, Y).
#' @param proteinName, optional character string. NA by default.
#'   Used to either add the name of the protein to the plot title.
#' @param ... additional arguments to be passed to
#'   \code{\link{structuralTendency}} and
#'   \code{\link[ggplot2]{ggplot}}
#' @param disorderPromoting,disorderNeutral,orderPromoting character vectors
#'    of individual residues to be matched with the input sequence. Defaults:
#'    \itemize{
#'      \item disorderPromoting = c("P", "E", "S", "Q", "K", "A", "G")
#'      \item orderPromoting =
#'         c("M", "N", "V", "H", "L", "F", "Y", "I", "W", "C")
#'      \item disorderNeutral = c("D", "T", "R")
#'    }
#'    It is not recommended to change these.
#' @return a data frame containing each residue from the sequence
#'   matched with its structural tendency, defined by disorderPromoting,
#'   disorderNeutral, and orderPromoting.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @family structural tendency
#' @references
#'   Uversky, V. N. (2013). A decade and a half of protein intrinsic disorder:
#'   Biology still waits for physics. Protein Science, 22(6), 693-724.
#'   \url{https://doi.org/10.1002/pro.2261}. \cr
#'   Kulkarni, Prakash, and Vladimir N. Uversky. "Intrinsically
#'   disordered proteins: the dark horse of the dark proteome."
#'   Proteomics 18.21-22 (2018): 1800061.
#'   \url{https://doi.org/10.1002/pmic.201800061}.
#' @export
#' @examples
#' #Amino acid sequences can be character strings
#' aaString <- "ACDEFGHIKLMNPQRSTVWY"
#' #Amino acid sequences can also be character vectors
#' aaVector <- c("A", "C", "D", "E", "F",
#'               "G", "H", "I", "K", "L",
#'               "M", "N", "P", "Q", "R",
#'               "S", "T", "V", "W", "Y")
#' #Alternativly, .fasta files can also be used by providing
#' ##The path to the file as a character string

#' structuralTendencyPlot(aaString)
#' structuralTendencyPlot(aaVector)
#'
#' #The plot can be a pie chart (default)
#' structuralTendencyPlot(aaString,
#'                     graphType = "pie")
#'
#' #Or the plot can be a bar graph
#' structuralTendencyPlot(aaString,
#'                     graphType = "bar")
#'
#' #To display general tendency rather than residues, set summarize = T
#' structuralTendencyPlot(aaString,
#'                     graphType = "pie",
#'                     summarize = TRUE)
#'
#' structuralTendencyPlot(aaString,
#'                     graphType = "bar",
#'                     summarize = TRUE)
#'
#' #If you wish to export this as a dataframe, set graphType = "none"
#' exampleDF <- structuralTendencyPlot(aaString,
#'                                   graphType = "none")
#' head(exampleDF)
#'
#' #If using a different definition of disordered residues
#' ##These residues are labeled as such from Dunker et al (2001),
#' ##"Intrinsically disordered protein."
#' structuralTendencyPlot(aaString,
#'               disorderPromoting = c("A", "R", "G", "Q", "S", "P", "E", "K"),
#'               disorderNeutral = c("H", "M", "T", "D"),
#'               orderPromoting = c("W", "C", "F", "I", "Y", "V", "L", "N"),
#'               graphType = "bar",
#'               alphabetical = TRUE)

structuralTendencyPlot <- function(
  sequence,
  graphType = "pie",
  summarize = FALSE,
  proteinName = NA,
  alphabetical = FALSE,
  disorderPromoting = c("P", "E", "S", "Q", "K", "A", "G"),
  disorderNeutral = c("D", "T", "R"),
  orderPromoting = c("M", "N", "V", "H", "L", "F", "Y", "I", "W", "C"),
  ...) {  #How to order output results (summarize = F)

  #---

  if (is.logical(summarize) == FALSE ||
      is.logical(alphabetical) == FALSE) {
    stop("summarize and alphabetical must be logical values")
  }
  if (!graphType %in% c("pie", "bar", "none")) {
    stop('invalid argument for graphType.
       Set graphType = c("pie", "bar", "none")')
  }

  structuralTendencyDF <- structuralTendency(sequence = sequence,
                                          disorderPromoting = disorderPromoting,
                                          disorderNeutral = disorderNeutral,
                                          orderPromoting = orderPromoting)
  sequenceLength <- nrow(structuralTendencyDF)

  if (summarize) {
    structuralTendencyDF <- data.frame(table(structuralTendencyDF$Tendency))
    names(structuralTendencyDF) <- c("Tendency", "Total")
    structuralTendencyDF$Frequency <- structuralTendencyDF$Total /
      sequenceLength * 100
    structuralTendencyDF$AA <- as.character(structuralTendencyDF$Tendency)

  } else {

    structuralTendencyDF <- structuralTendencyDF[, 2:3]
    residueFrequencyDF <- data.frame(table(structuralTendencyDF$AA))
    names(residueFrequencyDF) <- c("AA", "Total")
    structuralTendencyDF <- unique(structuralTendencyDF)
    structuralTendencyDF <- merge(structuralTendencyDF, residueFrequencyDF)
    structuralTendencyDF$Frequency <- round(structuralTendencyDF$Total /
                                              sequenceLength * 100, 3)
    names(structuralTendencyDF) <- c("AA", "Tendency", "Total", "Frequency")

    if (!alphabetical) {
      aaOrder <- c("P", "E", "S", "Q", "K", "A", "G",
                   "D", "T", "R",
                   "M", "N", "V", "H", "L", "F", "Y", "I", "W", "C")
      structuralTendencyDF$AA <- factor(structuralTendencyDF$AA,
                                        levels = aaOrder)
    }
  }
  if (!graphType == "none") {
    if (graphType == "bar") {
      gg <- ggplot2::ggplot(data = structuralTendencyDF,
                            ggplot2::aes_(x = ~ AA,
                                          y = ~ Frequency,
                                          fill = ~ Tendency,
                                          group = ~ Tendency)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_bw()
    }
    if (graphType == "pie") {
      #------Data needs mutated to label residues
      structuralTendencyDF <- structuralTendencyDF %>%
        dplyr::arrange(dplyr::desc(.data$Tendency)) %>%
        dplyr::mutate(prop = .data$Total / sum(.data$Total) * 100) %>%
        dplyr::mutate(ypos = cumsum(.data$prop) - 0.5 * .data$prop)
      gg <- ggplot2::ggplot(structuralTendencyDF,
                            ggplot2::aes_(x = "",
                                          y = ~ prop,
                                          fill = ~ Tendency)) +
        ggplot2::geom_bar(stat = "identity",
                          width = 1,
                          color = "white") +
        ggplot2::coord_polar("y",
                             start = 0) +
        ggplot2::theme_void() +
        ggplot2::geom_text(ggplot2::aes_(y = ~ ypos,
                                         label = ~ AA),
                           color = "white",
                           size = 4)
    }

    if (!is.na(proteinName)) {
      plotTitle <- paste("Compositional Profile of ",
                         proteinName,
                         sep = "", collapse = "")
    } else {
      plotTitle <- "Compositional Profile"
    }

    yTitle <- "Amino Acid Composition (as % of sequence length)"
    gg <- gg + ggplot2::labs(title = plotTitle,
                             y = yTitle)
    gg <- gg +
      ggplot2::theme(legend.position = "top",
                     plot.title = ggplot2::element_text(hjust = 0.5))
    return(gg)
  } else {
    return(structuralTendencyDF)
  }
}
