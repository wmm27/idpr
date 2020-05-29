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
#' @importFrom magrittr %>%
#' @importFrom rlang .data
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
        dplyr::arrange(desc(.data$Tendency)) %>%
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
