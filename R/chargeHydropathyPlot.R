#' Charge-Hydropathy Plot
#'
#' This function calculates the average net charge <R> and the average
#'   scaled hydropathy <H> and visualizes the data. There are known boundaries
#'   on the C-H plot that separate extended and collapsed proteins. \cr
#'   This was originally described in Uversky et al. (2000)\cr
#'   \url{https://doi.org/10.1002/1097-0134(20001115)41:3<415::AID-PROT130>3.0.CO;2-7}
#'   . \cr
#'   The plot returned is based on the charge-hydropathy plot from
#'   Uversky (2016) \url{https://doi.org/10.1080/21690707.2015.1135015}. \cr
#'   See Uversky (2019) \url{https://doi.org/10.3389/fphy.2019.00010} for
#'   additional information and a recent review on the topic.
#'   This plot has also been referred to as a "Uversky Plot".
#' @param sequence amino acid sequence (or pathway to a fasta file)
#'   as a character string. Supports multiple sequences / files, as a
#'   character vector of strings. Additionally, this supports a single protein
#'   as character vectors. Multiple proteins are not supported as a character
#'   vector of single characters.
#' @param displayInsolubility logical value, TRUE by default.
#'   This adds (or removes when FALSE) the vertical line
#'   separating collapsed proteins and insoluble proteins
#' @param insolubleValue numerical value. 0.7 by default.
#'   Ignored when \code{displayInsolubility = FALSE}. Plots the vertical line
#'   \eqn{<H> = displayInsolubility}.
#' @param proteinName,customPlotTitle optional character string. NA by default.
#'   Used to either add the name of the protein to the plot title when there
#'   is only one protein, or to create a custom plot title for the output.
#' @param pKaSet pKa set used for charge calculations. See
#'   \code{\link{netCharge}} for additional details
#' @param pH numeric value, 7.0 by default.
#'   The environmental pH used to calculate residue charge.
#' @param ... additional arguments to be passed to
#'   \link[idpr:netCharge]{idpr::netCharge()},
#'   \link[idpr:meanScaledHydropathy]{idpr::meanScaledHydropathy()} or
#'   \code{\link[ggplot2]{ggplot}}
#' @importFrom ggplot2 aes aes_
#' @return Graphical values of Charge-Hydropathy Plot
#' @seealso \code{\link{netCharge}} and
#'   \code{\link{meanScaledHydropathy}}
#'   for functions used to calculate values.
#' @references
#'   Kozlowski, L. P. (2016). IPC – Isoelectric Point Calculator. Biology
#'   Direct, 11(1), 55. \url{https://doi.org/10.1186/s13062-016-0159-9} \cr
#'   Kyte, J., & Doolittle, R. F. (1982). A simple method for
#'   displaying the hydropathic character of a protein.
#'   Journal of molecular biology, 157(1), 105-132. \cr
#'   Uversky, V. N. (2019). Intrinsically Disordered Proteins and Their
#'   “Mysterious” (Meta)Physics. Frontiers in Physics, 7(10).
#'   \url{https://doi.org/10.3389/fphy.2019.00010} \cr
#'   Uversky, V. N. (2016). Paradoxes and wonders of intrinsic disorder:
#'   Complexity of simplicity. Intrinsically Disordered Proteins, 4(1),
#'   e1135015. \url{https://doi.org/10.1080/21690707.2015.1135015} \cr
#'   Uversky, V. N., Gillespie, J. R., & Fink, A. L. (2000).
#'   Why are “natively unfolded” proteins unstructured under physiologic
#'   conditions?. Proteins: structure, function, and bioinformatics, 41(3),
#'   415-427.
#'   \url{https://doi.org/10.1002/1097-0134(20001115)41:3<415::AID-PROT130>3.0.CO;2-7}
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
#' chargeHydropathyPlot(sequence = aaString)
#' chargeHydropathyPlot( sequence = aaVector)
#'
#' #This function also supports multiple sequences
#' #only as character strings or .fasta files
#' multipleSeq <- c("ACDEFGHIKLMNPQRSTVWY",
#'                "ACDEFGHIK",
#'                "LMNPQRSTVW")
#' chargeHydropathyPlot(sequence = multipleSeq)
#'
#' #since it is a ggplot, we can add additional annotations or themes
#' chargeHydropathyPlot(
#'  sequence = multipleSeq)  +
#'   ggplot2::theme_void()
#'
#' chargeHydropathyPlot(
#'   sequence = multipleSeq)  +
#'   ggplot2::geom_hline(yintercept = 0,
#'                      color = "red")
#'
#' #choosing the pKa set used for calculations
#' chargeHydropathyPlot(
#'   sequence = multipleSeq,
#'   pKaSet = "EMBOSS")
#'


chargeHydropathyPlot <- function(
    sequence,
    displayInsolubility = TRUE,
    insolubleValue = 0.7,
    proteinName = NA,
    customPlotTitle = NA,
    pH = 7.0,
    pKaSet = "IPC_protein",
    ...) {

    if (nchar(sequence[1]) == 1) {
        sequence <- paste(sequence, sep = "", collapse = "")
    }
    #--- Calculating the C-H data for each protein
    nSequences <- length(sequence)
    dataCollected <- data.frame(matrix(nrow = nSequences,
                                        ncol = 3))
    names(dataCollected) <- c("sequence",
                            "avg_scaled_hydropathy",
                            "avg_net_charge")

    for (i in seq_len(nSequences)) {
        sequence.i <- sequence[i]
        dataCollected$sequence[i] <- sequence.i
        dataCollected$avg_scaled_hydropathy[i] <-
            meanScaledHydropathy(sequence = sequence.i)
        dataCollected$avg_net_charge[i] <-
            netCharge(sequence = sequence.i,
                        pKaSet = pKaSet,
                        pH = pH,
                        includeTermini = TRUE,
                        averaged = TRUE)
    }

    # ---- Math for plotting lines
    #The equations for the lines are:
    #  Boundary seperating IDPs and compact proteins
    #   <R> = 2.785 * <H> - 1.151
    #   <R> = -2.785 * <H> + 1.151
    #  Limits of CH space
    #   <R> = 1.125 * <H> - 1.125
    #   <R> = 1.000 - <H>
    #  Insoluble line
    #   <H> = 0.700 (or custom value)

    intersectionPointX <- (1.151 * 2) / (2.785 * 2)

    positiveBoundaryX <- (1.151 + 1) / (1 + 2.785)
    positiveBoundaryY <- (-1 * positiveBoundaryX) + 1

    negativeBoundaryX <- (1.151 + 1.125) / (1.125 + 2.785)
    negativeBoundaryY <- 1.125 * negativeBoundaryX - 1.125

    # --- making the ggplot
    gg <- ggplot2::ggplot(dataCollected,
                        aes_(x = ~ avg_scaled_hydropathy, y = ~ avg_net_charge))
    gg <- gg + ggplot2::geom_segment(aes(x = intersectionPointX,
                                        y = 0,
                                        xend = positiveBoundaryX,
                                        yend = positiveBoundaryY))
    gg <- gg + ggplot2::geom_segment(aes(x = intersectionPointX,
                                        y = 0,
                                        xend = negativeBoundaryX,
                                        yend = negativeBoundaryY))
    if (displayInsolubility) {
        if (!is.numeric(insolubleValue)) {
            stop("insolubleValue must be a numeric value.")
        }
        insolubleMax <- (-1 * insolubleValue) + 1
        insolubleMin <-  (1.125 * insolubleValue) - 1.125
        gg <- gg + ggplot2::geom_segment(aes(x = insolubleValue,
                                            y = insolubleMax,
                                            xend = insolubleValue,
                                            yend = insolubleMin))
    gg <- gg +
        ggplot2::geom_label(aes(x = 0.85, y = 0.35,
                                label = "Insoluble Proteins")) +
        ggplot2:: geom_label(aes(x = 0.7, y = 0.5,
                                label = "Collapsed Proteins"))
    } else {
        gg <- gg +
            ggplot2::geom_label(ggplot2::aes(x = 0.8,
                                            y = 0.4,
                                            label = "Collapsed Proteins"))
    }

    gg <- gg +
        ggplot2::geom_label(ggplot2::aes(x = 0.4,
                                        y = 0.8,
                                        label = "Extended IDPs"))

    #Values cannot exceede the logical space within the C-H plot
    gg <- gg +
        ggplot2::geom_segment(ggplot2::aes(x = 1, y = 0,
                                            xend = 0, yend = 1)) +
        ggplot2::geom_segment(ggplot2::aes(x = 1, y = 0,
                                            xend = 0, yend = -1.125)) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, y = 1,
                                            xend = 0, yend = -1.125))
    xLabel <- paste("Mean Scaled Hydropathy")
    yLabel <- paste("Mean Net Charge")

    if (is.na(customPlotTitle)) {
        if (nSequences == 1 &&
            !is.na(proteinName)) {
            ggTitle <- paste("Charge-Hydropathy Plot of ", proteinName,
                            sep = "", collapse = "")
        }
        if (nSequences > 1 ||
            is.na(proteinName)) {
            ggTitle <- "Charge-Hydropathy Plot"
        }
    } else {
        ggTitle <- customPlotTitle
    }
    gg <- gg + ggplot2::geom_point(color = "#92140C") +
        ggplot2::theme_minimal() + ggplot2::xlim(0, 1) +
        ggplot2::ylim(-1.125, 1) + ggplot2::labs(y = yLabel, x = xLabel,
                                                    title = ggTitle)
    return(gg)
}
