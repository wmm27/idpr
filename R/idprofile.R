#' IDPRofile From idpr Package
#'
#' The IDPRofile is a summation of many features of the idpr package,
#'   conveniently grouped into one function for quick analysis. This combines
#'   many plotting functions in this package. These include:\cr
#'   \code{\link{chargeHydropathyPlot}}\cr
#'   \code{\link{chargeCalculationLocal}}\cr
#'   \code{\link{scaledHydropathyLocal}}\cr
#'   \code{\link{structuralTendencyPlot}}\cr
#'   All of the above linked functions only require the sequence argument
#'   to ouput plots of characteristics associated with IDPs. The function also
#'   includes options for IUPred functions. The function does one of the
#'   following based on user-specified parameters:\cr
#'   \code{\link{iupred}}\cr
#'   \code{\link{iupredAnchor}}\cr
#'   \code{\link{iupredRedox}}\cr
#'   The IUPred function used depends on the argument of iupredType. All
#'   require the UniProt Accession to make a proper connection to the IUPred2A
#'   REST API. If the UniProt Accession is not specified, the IUPred plot is
#'   skipped.
#' @param sequence amino acid sequence as a single character string or vector of
#'    single characters. It also supports a single character string that
#'    specifies the locaion of a .fasta or .fa file.
#' @param uniprotAccession character string specifying the UniProt Accession of
#'   the protein of interest. Used to fetch predictions from IUPreds REST API.
#'   Default is NA. Keep as NA if you do not have a UniProt Accession.
#'
#' @param proteinName character string, optional.
#'   Used to add protein name to the title in ggplot.
#' @inheritParams chargeCalculationLocal
#' @inheritParams structuralTendencyPlot
#' @param structuralTendencyType a character string specifying the type of plot
#'   the \code{\link{structuralTendencyPlot}} should output.
#'   Can be "bar" or "pie".
#'   Equivilent argument to graphType= in the linked function.
#'   "bar" by default.
#' @param structuralTendencySummarize a logical value specifying the
#'    \code{\link{structuralTendencyPlot}} should be summarized into broad
#'    categories. Equivilent argument to summarize= in the linked function.
#'    FALSE by defualt
##' @param disorderPromoting,disorderNeutral,orderPromoting character vectors
#'    of individual residues to be matched with the input sequence. Defaults:
#'    \itemize{
#'      \item disorderPromoting = c("P", "E", "S", "Q", "K", "A", "G")
#'      \item orderPromoting =
#'         c("M", "N", "V", "H", "L", "F", "Y", "I", "W", "C")
#'      \item disorderNeutral = c("D", "T", "R")
#'    }
#'    It is not reccomended to change these. Arguments passed to
#'    \code{\link{structuralTendencyPlot}}
#' @param iupredType character string specifying the type of IUPred2 prediction
#'   to retrieve. Can be c("long", "short", "glob", "anchor", "redox"). "long"
#'   by default. "long", "short", and "glob" use the \code{\link{iupred}}
#'   function and specify the type of plot. Both "redox" and "anchor" use "long"
#'   for predictions, but are context dependent. "anchor" uses
#'   \code{\link{iupredAnchor}} to get predictions of disorder with IUPred2 and
#'   preditions of induced folding based on ANCHOR2 predictions (Shown with a
#'   red line). "redox" uses \code{\link{iupredRedox}} to make predictions of
#'   disorder based on enviornmental conditions. Regions of predicted
#'   environmental sensitivity are highlighted. See the respetive functions
#'   for more details. This is skipped if uniprotAccession = NA.
#' @return 4 or 5 plots, depending if a UniProt Accession is provided.
#' @export
#' @seealso \code{\link{chargeHydropathyPlot}}\cr
#'   \code{\link{chargeCalculationLocal}}\cr
#'   \code{\link{scaledHydropathyLocal}}\cr
#'   \code{\link{structuralTendencyPlot}}\cr
#'   \code{\link{iupred}}\cr
#'   \code{\link{iupredAnchor}}\cr
#'   \code{\link{iupredRedox}}\cr
#'
#' @section Citations for each Plot:
#'   \itemize{
#'     \item \code{\link{chargeHydropathyPlot}}
#'       \itemize{
#'         \item{Kozlowski, L. P. (2016). IPC – Isoelectric Point Calculator.
#'               Biology Direct, 11(1), 55.
#'               https://doi.org/10.1186/s13062-016-0159-9}
#'         \item{Kyte, J., & Doolittle, R. F. (1982).
#'               A simple method for displaying the hydropathic character
#'               of a protein. Journal of molecular biology, 157(1), 105-132.}
#'         \item{Uversky, V. N. (2019). Intrinsically Disordered Proteins and
#'               Their “Mysterious” (Meta)Physics. Frontiers in Physics, 7(10).
#'               https://doi.org/10.3389/fphy.2019.00010}
#'         \item{If a pKa set is specified, see \code{\link{pKaData}}}
#'     }
#'     \item \code{\link{chargeCalculationLocal}}
#'       \itemize{
#'         \item{Kozlowski, L. P. (2016). IPC – Isoelectric Point Calculator.
#'               Biology Direct, 11(1), 55.
#'               https://doi.org/10.1186/s13062-016-0159-9}
#'         \item{If a pKa set is specified, see \code{\link{pKaData}}}
#'     }
#'     \item \code{\link{scaledHydropathyLocal}}
#'       \itemize{
#'         \item{Kyte, J., & Doolittle, R. F. (1982).
#'               A simple method for displaying the hydropathic character
#'               of a protein. Journal of molecular biology, 157(1), 105-132.}
#'     }
#'     \item \code{\link{structuralTendencyPlot}}
#'       \itemize{
#'         \item{Kulkarni, Prakash, and Vladimir N. Uversky.
#'         "Intrinsically disordered proteins: the dark horse of the
#'          dark proteome." Proteomics 18.21-22 (2018): 1800061.
#'          https://doi.org/10.1002/pmic.201800061.}
#'     }
#'     \item \code{\link{iupred}},
#'           \code{\link{iupredAnchor}},
#'           \code{\link{iupredRedox}}
#'       \itemize{
#'         \item{Bálint Mészáros, Gábor Erdős, Zsuzsanna Dosztányi, IUPred2A:
#'         context-dependent prediction of protein disorder as a function of
#'         redox state and protein binding, Nucleic Acids Research, Volume 46,
#'          Issue W1, 2 July 2018, Pages W329–W337,
#'          https://doi.org/10.1093/nar/gky384}
#'         \item{Erdős, G., & Dosztányi, Z. (2020). Analyzing protein disorder
#'         with IUPred2A. Current Protocols in Bioinformatics, 70, e99.
#'         https://doi.org/10.1002/cpbi.99}
#'     }
#'  }
#' @examples
#' \dontrun{
#' #For most functions, a protein sequence is all that is needed.
#'
#' #The UniProt ID is optional but reccomened for IUPred results.
#' proteinID <- "P04637"
#' iupredDF <- iupred(proteinID,
#'                    plotResults = FALSE)
#' iupredSeq <- iupredDF$AA
#'
#' idprofile(
#'   sequence = iupredSeq,
#'   uniprotAccession = proteinID)
#'
#'
#' #changing the iupred to redox
#' ## and getting a pie chart for structuralTendency.
#' idprofile(
#'   sequence = iupredSeq,
#'   uniprotAccession = proteinID,
#'   iupredType = "redox",
#'   structuralTendencyType = "pie")
#' }

idprofile <- function(sequence,
                      uniprotAccession = NA,
                      proteinName = NA,
                      window = 9,
                      pH = 7.2,
                      pKaSet = "IPC_protein",
                      structuralTendencyType = "bar",
                      structuralTendencySummarize = FALSE,
                      disorderPromoting = c("P", "E", "S", "Q", "K", "A", "G"),
                      disorderNeutral = c("D", "T", "R"),
                      orderPromoting = c("M", "N", "V", "H", "L",
                                         "F", "Y", "I", "W", "C"),
                      iupredType = "long"
) {

  #---- General Plots
  rhPlot <- chargeHydropathyPlot(sequence = sequence,
                                 pH = pH,
                                 pKaSet = pKaSet,
                                 proteinName = proteinName)
  chargePlot <- chargeCalculationLocal(sequence = sequence,
                                       window = window,
                                       plotResults = TRUE,
                                       pH = pH,
                                       proteinName = proteinName)
  hydropPlot <- scaledHydropathyLocal(sequence = sequence,
                                      window = window,
                                      plotResults = TRUE,
                                      pKaSet = pKaSet,
                                      proteinName = proteinName)
  tendencyPlot <- structuralTendencyPlot(sequence = sequence,
                                        graphType = structuralTendencyType,
                                        summarize = structuralTendencySummarize,
                                        disorderPromoting = disorderPromoting,
                                        disorderNeutral = disorderNeutral,
                                        orderPromoting = orderPromoting,
                                        proteinName = proteinName)
  #-------- Adding IUPred Plot based on which type
  if (!is.na(uniprotAccession)) {
    if (iupredType %in% c("long", "short", "glob")) {
      iupredPlot <- iupred(uniprotAccession,
                           iupredType = iupredType,
                           plotResults = TRUE,
                           proteinName = proteinName)
    }
    if (iupredType == "anchor") {
      iupredPlot <- iupredAnchor(uniprotAccession,
                                 plotResults = TRUE,
                                 proteinName = proteinName)
    }
    if (iupredType == "redox") {
      iupredPlot <- iupredRedox(uniprotAccession,
                                plotResults = TRUE,
                                proteinName = proteinName)
    }
  } else {
    iupredPlot <- ggplot2::ggplot() +
      ggplot2::annotate("text",
                        x = 1,
                        y = 1,
                        label = "No Uniprot Accession provided.\n
                                   IUPred plot skipped") +
      ggplot2::theme_void()
  }
  plotList <- list((rhPlot),
                   (tendencyPlot),
                   (chargePlot),
                   (hydropPlot),
                   (iupredPlot))
  return(plotList)
}
