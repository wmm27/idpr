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
#'   When \code{graphType = "none"}, the output is the dataframe that would
#'   otherwise be used to plot the data.
#' @param summarize logical value, FALSE by default.
#'   When \code{summarize = TRUE}, each residue is agregated into Disorder
#'   Tendancy Groups. (See \code{\link{structuralTendency}} for more details).
#'   When \code{summarize = FALSE}, residue identity is preserved and
#'   the output is colored by Disorder Tendancy Groups.
#' @param alphabetical logical value, FALSE by default.
#'   Order of residues on plot axis. Only relevant when
#'    \code{summarize = FALSE}, otherwise is ignored.
#'    If FALSE, ordering is grouped by Disorder Tendancy (P, E, S, ..., W, C).
#'    If TRUE, the residues are ordered alphabetically (A, C, D, E, ..., W, Y).
#' @param proteinName, optional character string. NA by default.
#'   Used to either add the name of the protein to the plot title.
#' @param ... additional arguments to be passed to
#'   \code{\link{structuralTendency}} and
#'   \code{\link[ggplot2]{ggplot}}
#'
#' @return a dataframe containing each residue from the sequence
#'   matched with its structural tendancy, defined by disorderPromoting,
#'   disorderNeutral, and orderPromoting.
#' @importFrom magrittr %>%
#' @family structural tendency
#' @references Kulkarni, Prakash, and Vladimir N. Uversky. "Intrinsically
#'   disordered proteins: the dark horse of the dark proteome."
#'   Proteomics 18.21-22 (2018): 1800061.
#'   \url{https://doi.org/10.1002/pmic.201800061}.
#' @export

structuralTendencyPlot <- function(
  sequence,
  graphType = "pie",
  summarize = FALSE,
  proteinName = NA,
  alphabetical = FALSE,
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


  if (summarize) {
    structuralTendencyDF <- structuralTendency(sequence = sequence)
    sequenceLength <- nrow(structuralTendencyDF)
    structuralTendencyDF <- data.frame(table(structuralTendencyDF$Tendency))
    names(structuralTendencyDF) <- c("Tendency", "Total")
    structuralTendencyDF$Frequency <- structuralTendencyDF$Total /
      sequenceLength * 100
    structuralTendencyDF$AA <- as.character(structuralTendencyDF$Tendency)

  } else {
    seqCharacterVector <- sequenceCheck(
      sequence = sequence,
      sequenceName = NA,
      method = "Stop",
      outputType = "Vector",
      supressOutputMessage = T)
    sequenceLength <- length(seqCharacterVector)

    residueFrequencyDF <- data.frame(table(seqCharacterVector))
    names(residueFrequencyDF) <- c("AA", "Total")
    residueVector <- as.character(residueFrequencyDF$AA)
    structuralTendencyDF <- structuralTendency(sequence = residueVector)
    removeVector <- names(structuralTendencyDF) == "Position"
    structuralTendencyDF <- structuralTendencyDF[, !removeVector]
    structuralTendencyDF <- merge(structuralTendencyDF, residueFrequencyDF)
    structuralTendencyDF$Frequency <- round(structuralTendencyDF$Total /
                                              sequenceLength * 100, 3)

    if (alphabetical == F) {
      aaOrder <- c("P", "E", "S", "Q", "K", "A", "G",
                   "D", "T", "R",
                   "M", "N", "V", "H", "L", "F", "Y", "I", "W", "C")
      structuralTendencyDF$AA <- factor(structuralTendencyDF$AA,
                                        levels = aaOrder)
    }
  }
  if (!graphType == "none") {
    if (graphType == "bar") {
      gg <- ggplot2::ggplot(structuralTendencyDF,
                            ggplot2::aes(x = AA,
                                         y = Frequency,
                                         fill = Tendency,
                                         group = Tendency))

      gg <- gg + ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_bw()
    }
    if (graphType == "pie") {
      #------Data needs mutated to label residues
      structuralTendencyDF <- structuralTendencyDF %>%
        dplyr::arrange(desc(Tendency)) %>%
        dplyr::mutate(prop = Total / sum(structuralTendencyDF$Total) * 100) %>%
        dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop)


      gg <- ggplot2::ggplot(structuralTendencyDF,
                            ggplot2::aes(x = "",
                                         y = prop,
                                         fill = Tendency))

      gg <- gg + ggplot2::geom_bar(stat = "identity",
                                   width = 1,
                                   color = "white") +
        ggplot2::coord_polar("y",
                             start = 0) +
        ggplot2::theme_void() +
        ggplot2::geom_text(ggplot2::aes(y = ypos,
                                        label = AA),
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

    plot(gg)
  } else {
    return(structuralTendencyDF)
  }
}
