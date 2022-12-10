#' Prediction of Intrinsic Disorder with FoldIndex method in R
#'
#' This is used to calculate the prediction of intrinsic disorder based on
#'   the scaled hydropathy and absolute net charge of an amino acid
#'   sequence using a sliding window. FoldIndex described this relationship and
#'   implemented it graphically in 2005 by Prilusky, Felder, et al,
#'   and this tool has been implemented
#'   into multiple disorder prediction programs. When windows have a negative
#'   score (<0) sequences are predicted as disordered.
#'   When windows have a positive score (>0) sequences are predicted as
#'   disordered. Graphically, this cutoff is displayed by the dashed
#'   line at y = 0. Calculations are at pH 7.0 based on the described method and
#'   the default is a sliding window of size 51.
#'
#'   The output is either a data frame or graph
#'   showing the calculated scores for each window along the sequence.
#'   The equation used was originally described in Uversky et al. (2000)\cr
#'   \url{https://doi.org/10.1002/1097-0134(20001115)41:3<415::AID-PROT130>3.0.CO;2-7}
#'   . \cr
#'   The FoldIndex method of using a sliding window and utilizing the Uversky
#'   equation is described in Prilusky, J., Felder, C. E., et al. (2005). \cr
#'   FoldIndex: a simple tool to predict whether a given protein sequence \cr
#'   is intrinsically unfolded. Bioinformatics, 21(16), 3435-3438. \cr
#'
#' @inheritParams sequenceCheck
#' @inheritParams chargeCalculationLocal
#' @param window a positive, odd integer. 51 by default.
#'   Sets the size of sliding window, must be an odd number.
#'   The window determines the number of residues to be analyzed and averaged
#'   for each position along the sequence.
#' @param plotResults logical value, TRUE by default.
#'   If \code{plotResults = TRUE} a plot will be the output.
#'   If \code{plotResults = FALSE} the output is a data frame with scores for
#'   each window analyzed.
#' @param proteinName character string with length = 1.
#'   optional setting to replace the name of the plot if plotResults = TRUE.
#' @param ... any additional parameters, especially those for plotting.
#' @return see plotResults argument
#' @family scaled hydropathy functions
#' @seealso \code{\link{KDNorm}} for residue hydropathy values.
#'   See \code{\link{pKaData}} for residue pKa values and citations. See
#'   \code{\link{hendersonHasselbalch}} for charge calculations.

#' @section Plot Colors:
#'   For users who wish to keep a common aesthetic, the following colors are
#'   used when plotResults = TRUE. \cr
#'   \itemize{
#'   \item Dynamic line colors: \itemize{
#'   \item Close to -1 = "#9672E6"
#'   \item Close to 1 = "#D1A63F"
#'   \item Close to midpoint = "grey65" or "#A6A6A6"}}
#'
#' @references
#'   Kyte, J., & Doolittle, R. F. (1982). A simple method for
#'   displaying the hydropathic character of a protein.
#'   Journal of molecular biology, 157(1), 105-132.
#'   Kozlowski, L. P. (2016). IPC – Isoelectric Point Calculator. Biology
#'   Direct, 11(1), 55. \url{https://doi.org/10.1186/s13062-016-0159-9} \cr
#'   Kyte, J., & Doolittle, R. F. (1982). A simple method for
#'   displaying the hydropathic character of a protein.
#'   Journal of molecular biology, 157(1), 105-132. \cr
#'   Prilusky, J., Felder, C. E., et al. (2005). \cr
#'   FoldIndex: a simple tool to predict whether a given protein sequence \cr
#'   is intrinsically unfolded. Bioinformatics, 21(16), 3435-3438. \cr
#'   Uversky, V. N., Gillespie, J. R., & Fink, A. L. (2000).
#'   Why are “natively unfolded” proteins unstructured under physiologic
#'   conditions?. Proteins: structure, function, and bioinformatics, 41(3),
#'   415-427.
#'   \url{https://doi.org/10.1002/1097-0134(20001115)41:3<415::AID-PROT130>3.0.CO;2-7}
#' @export

foldIndexR <- function(sequence,
                       window = 51,
                       proteinName = NA,
                       pKaSet = "IPC_protein",
                       plotResults = TRUE,
                       ...) {
    chargeDF <-
        chargeCalculationLocal(sequence = sequence, window = window,
                               pH = 7.0, pKaSet = pKaSet,
                               plotResults = FALSE)
    chargeDF$scaledWindowCharge <- chargeDF$windowCharge / window
    hydropDF <-  scaledHydropathyLocal(sequence = sequence,
                                       window = window,
                                       plotResults = FALSE)
    mergeDF <- merge(hydropDF, chargeDF)
    mergeDF$foldIndex <-
        mergeDF$WindowHydropathy * 2.785 -
        abs(mergeDF$scaledWindowCharge) - 1.151
    if (plotResults) {
        plotTitle <- "FoldIndex Prediction of Intrinsic Disorder"
        if (!is.na(proteinName)) {
            plotTitle <-
                paste0("FoldIndex Prediction of Intrinsic Disorder in ",
                       proteinName, sep = "")
        }
        gg <-  sequencePlot(position = mergeDF$Position,
                            property = mergeDF$foldIndex,
                            hline = 0, dynamicColor = mergeDF$foldIndex,
                            customColors = c("#9672E6", "#D1A63F", "grey65"),
                            customTitle = NA, propertyLimits = c(-1, 1))
        gg <- gg + ggplot2::labs(title = plotTitle, y = "Score")
        return(gg)
    } else {
        return(mergeDF)
    }
}
