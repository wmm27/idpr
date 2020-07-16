#' Prediction of Intrinsic Disorder with IUPred2A
#'
#' This function makes a connection to the IUPred2A REST API based on the type
#'   of analysis and UniProt accession number. This requires the user to know
#'   the accession number of their protein and a connection to the internet.
#'   The results are then formatted to match output in the idpr package. \cr \cr
#'   Predictions are made on a scale of 0-1, where any residues with a score
#'   over 0.5 are predicted to be disordered, and any residue scoring below 0.5
#'   are predicted to be ordered (when using "long" and "short" predictions).\cr
#'   The output is either a graph (ggplot) or data frame of predictions.
#'   \cr\cr
#'   \strong{iupred()} is used for standard predictions of intrinsic disorder
#'   of an amino acid sequence. This is the core of predictions.
#'   Predictions vary by iupredType (details below)
#'   The results are either a ggplot or data frame of the fetched IUPred2.
#'   predictions.
#'   \cr
#'   \strong{iupredAnchor()} is used to combine the output of IUPred2 long with
#'   ANCHOR2 predictions. ANCHOR2 is a context-dependent predictor of binding
#'   regions for protein-protein interactions. The results are either a ggplot
#'   with 2 lines, one for IUPred2 long and another for ANCHOR predictions, or
#'   a data frame with both IUPred2 long and ANCHOR Predictions. Values are
#'   fetched by the IUPred2A REST API.
#'   \cr
#'   \strong{iupredRedox()} is used to predict redox-sensitive regions that may
#'   experience induced folding upon changing environments.
#'   This is a context-dependent predictor of disordered regions depending on
#'   a reducing (plus) or oxidizing (minus) environment. The results can be
#'   a ggplot with two IUPred2 long predictions, one for plus and another for
#'   minus environments, with redox sensitive regions shaded (if predicted).
#'   Alternatively, the results can be a data frame with both IUPred2 long plus
#'   and minus predictions as well as a column of logical values where a residue
#'   that is TRUE is predicted to be in a redox sensitive region. Values are
#'   fetched by the IUPred2A REST API.
#'   \cr \cr
#'   IUPred2 website is located at \url{https://iupred2a.elte.hu/}.
#'   For detailed information on using IUPred2A, please refer to
#'   \href{https://doi.org/10.1002/cpbi.99}{Erdős & Dosztány (2020)}
#'   Analyzing protein disorder with IUPred2A.
#'   Current Protocols in Bioinformatics, 70, e99.
#'   Additionally, please see
#'   \href{https://doi.org/10.1093/nar/gky384}{Mészáros et al (2019)}
#'   for further information, theory, and applications of IUPred2A.
#'   \cr \cr
#'   \strong{Please cite these articles if you use any iupred function.}
#'
#' @param uniprotAccession character string specifying the UniProt Accession
#'   of the protein of interest. Used to fetch predictions from IUPreds REST API
#' @param iupredType character string. "long" by default. accepted types are
#'   c("long", "short", "glob"). See "Prediction Type" information below.
#' @param proteinName character string, optional. Used to add protein name
#'   to the title in ggplot. Ignored if \code{plotResults = FALSE}.
#' @param plotResults logical value. TRUE by default.
#'   If \code{plotResults = TRUE}, a ggplot of IUPred predictions is returned
#'   If \code{plotResults = FALSE}, a data frame of predictions is returned.
#' @return see plotResults argument.
#'
#' @section Prediction Type:
#'   Information from \url{https://iupred2a.elte.hu/help_new} on 5.22.20
#'   Additionally, see the sources for further details and source information.
#'   This is only relevant for iupred(). iupredAnchor() and iupredRedox()
#'   always utilize "long" for data in the REST API.
#'   \itemize{
#'     \item Long predictions of disorder (Default)
#'        \itemize{
#'          \item when iupredType = "long"
#'          \item Optimized for global predictions of disorder, specifically
#'            disordered regions over 30 amino acids in length.
#'          \item "long" is always used for iupredAnchor() and iupredRedox().
#'        }
#'      \item Short predictions of disorder
#'        \itemize{
#'          \item when iupredType = "short"
#'          \item Best for predicting small regions of disorder, especially
#'            in mostly structured proteins.
#'          \item Has adjustments for termini, since sequence ends are often
#'            disordered.
#'        }
#'      \item Structured predictions
#'        \itemize{
#'          \item when iupredType = "glob"
#'          \item Used to predict regions of globular folding.
#'          \item please see
#'            \href{https://doi.org/10.1002/cpbi.99}{Erdős & Dosztány (2020)}
#'            for further information on interpreting these results.
#'        }
#'    }
#' @source Bálint Mészáros, Gábor Erdős, Zsuzsanna Dosztányi,
#'   IUPred2A: context-dependent prediction of protein disorder as a function of
#'   redox state and protein binding, Nucleic Acids Research, Volume 46, Issue
#'   W1, 2 July 2018, Pages W329–W337, \url{https://doi.org/10.1093/nar/gky384}
#'   \cr\cr
#'   Erdős, G., & Dosztányi, Z. (2020). Analyzing protein disorder with
#'   IUPred2A. Current Protocols in Bioinformatics, 70, e99.
#'   \url{https://doi.org/10.1002/cpbi.99}
#' @export
#' @section Plot Colors:
#'   For users who wish to keep a common aesthetic, the following colors are
#'   used when plotResults = TRUE. \cr \itemize{
#'   \item iupred() iupredType = 'long', 'short', or 'glob'. Additionally,
#'      the 'long' prediction with iupredAnchor(). \itemize{
#'   \item Dynamic iupred line colors: \itemize{
#'   \item Close to 0 = "darkolivegreen3" or "#A2CD5A"
#'   \item Close to 1 = "darkorchid1" or "#BF3EFF"
#'   \item Close to 0.5 (midpoint) = "grey65" or "#A6A6A6"}}
#'   \item iupredAnchor : \itemize{
#'   \item Solid Line (ANCHOR2 Score) = "#92140C"}
#'   \item iupredRedox: \itemize{
#'   \item iupredPlus line = "darkorchid1" or "#BF3EFF"
#'   \item iupredMin line = "#348AA7"
#'   \item redox sensitive regions = "#5DD39E"}
#'   }
#' @name iupred
#' @examples
#' #A UniProt Accession must be specified.
#' ##this example uses human P53.
#' TP53_UniProt <- "P04637"
#' \dontrun{
#' #Getting data as a data frame
#' exampleDF_long <- iupred(uniprotAccession = TP53_UniProt,
#'                          iupredType = "long",
#'                          plotResults = FALSE)
#' head(exampleDF_long)
#'
#' exampleDF_short <- iupred(uniprotAccession = TP53_UniProt,
#'                           iupredType = "short",
#'                           plotResults = FALSE)
#' head(exampleDF_short)
#'
#' exampleDF_anchor <- iupredAnchor(uniprotAccession = TP53_UniProt,
#'                                  plotResults = FALSE)
#' head(exampleDF_anchor)
#'
#' exampleDF_redox <- iupredRedox(uniprotAccession = TP53_UniProt,
#'                               plotResults = FALSE)
#' head(exampleDF_redox)
#'
#' #Plotting
#'
#' iupred(uniprotAccession = TP53_UniProt,
#'     iupredType = "long",
#'     plotResults = TRUE)
#'
#' iupred(uniprotAccession = TP53_UniProt,
#'     iupredType = "short",
#'     plotResults = TRUE)
#'
#' iupredAnchor(uniprotAccession = TP53_UniProt,
#'           plotResults = TRUE)
#'
#' iupredRedox(uniprotAccession = TP53_UniProt,
#'          plotResults = TRUE)
#' }
#' #A valid internet connection is needed to make
#' ##A connection with the IUPred REST API

#----
iupred <- function(
    uniprotAccession,
    iupredType = "long",
    plotResults = TRUE,
    proteinName = NA) {
    #------
    #Connecting to IUPred2A REST API
    iupredURL <- paste("https://iupred2a.elte.hu/iupred2a/",
                        iupredType, "/", uniprotAccession, ".json", sep = "")
    iupredJson <- jsonlite::fromJSON(iupredURL)
    #-----
    #Reformatting data to be consistent in formatting across idpr
    iupredPrediction <- iupredJson$iupred2
    iupredSequence <- unlist(strsplit(iupredJson$sequence, ""))
    iupredSequence <- unlist(iupredSequence)
    seqLength <- length(iupredSequence)
    iupredDF <- data.frame(Position = seq_len(seqLength),
                            AA = iupredSequence,
                            IUPred2 = iupredPrediction)
    #Returning fetched results
    if (plotResults) {
        plotTitle <- "Prediction of Intrinsic Disorder"
        if (!is.na(proteinName)) {
            plotTitle <- paste("Prediction of Intrinsic Disorder in ",
                                proteinName,
                                sep = "")
        }
        plotSubtitle <- paste("By IUPred2A ", iupredJson$type, sep = "")
        gg <- sequencePlot(
                position = iupredDF$Position,
                property = iupredDF$IUPred2,
                hline = 0.5,
                dynamicColor = iupredDF$IUPred2,
                customColors = c("darkolivegreen3", "darkorchid1", "grey65"),
                customTitle = NA,
                propertyLimits = c(0, 1))

        gg <- gg + ggplot2::labs(title = plotTitle, subtitle = plotSubtitle)
        return(gg)
    } else {
        return(iupredDF)
    }
}

#' @rdname iupred
#' @export
iupredAnchor <- function(
    uniprotAccession,
    plotResults = TRUE,
    proteinName = NA) {
    #---- Connecting to IUPred2A REST API
    iupredURL <- paste("https://iupred2a.elte.hu/iupred2a/anchor/",
                        uniprotAccession, ".json", sep = "")
    iupredJson <- jsonlite::fromJSON(iupredURL)
    #----- Reformatting data to be consistent in formatting across idpr
    iupredPrediction <- iupredJson$iupred2
    anchorPrediction <- iupredJson$anchor2
    iupredSequence <- unlist(strsplit(iupredJson$sequence, ""))
    iupredSequence <- unlist(iupredSequence)
    seqLength <- length(iupredSequence)
    iupredDF <- data.frame(Position = seq_len(seqLength),
                            AA = iupredSequence,
                            IUPred2 = iupredPrediction,
                            ANCHOR2 = anchorPrediction)
    #------ Returning results
    if (plotResults) {
        plotTitle <- "Prediction of Intrinsic Disorder"
        if (!is.na(proteinName)) {
            plotTitle <- paste("Prediction of Intrinsic Disorder in ",
                                proteinName, sep = "")
        }
        plotSubtitle <- paste("By IUPred2A ", iupredJson$type,
                                " and ANCHOR2", sep = "")
        gg <- sequencePlot(
            position = iupredDF$Position,
            property = iupredDF$IUPred2,
            hline = 0.5,
            dynamicColor = iupredDF$IUPred2,
            customColors = c("darkolivegreen3", "darkorchid1", "grey65"),
            customTitle = NA,
            propertyLimits = c(0, 1))
        gg <-  gg + ggplot2::geom_line(data = iupredDF,
                    ggplot2::aes_(x = ~ Position,
                                y = ~ ANCHOR2),
                                color = "#92140C",
                                inherit.aes = FALSE)
        gg <- gg + ggplot2::labs(title = plotTitle, subtitle = plotSubtitle)
        return(gg)
    } else {
        return(iupredDF)
    }
}

#' @rdname iupred
#' @export
iupredRedox <-
    function(uniprotAccession, plotResults = TRUE, proteinName = NA) {
    iupredURL <- paste("https://iupred2a.elte.hu/iupred2a/redox/",
                        uniprotAccession, ".json", sep = "")
    iupredJson <- jsonlite::fromJSON(iupredURL)
    iupredPlus <- iupredJson$iupred2_redox_plus
    iupredMinus <- iupredJson$iupred2_redox_minus
    redoxSenstitiveDF <- as.data.frame(iupredJson$redox_sensitive_regions)
    iupredSequence <- unlist(strsplit(iupredJson$sequence, ""))
    seqLength <- length(iupredSequence)
    iupredDF <- data.frame(Position = seq_len(seqLength), AA = iupredSequence,
                            iupredPlus = iupredPlus, iupredMinus = iupredMinus)
    if (plotResults) {
        plotTitle <- "Prediction of Intrinsic Disorder"
        if (!is.na(proteinName)) {
            plotTitle <- paste("Prediction of Intrinsic Disorder in ",
                                proteinName, sep = "")
        }
        plotSubtitle <- paste("By IUPred2 ", iupredJson$type,
                            "|Based on Environmental Redox State", sep = "")
        gg <- ggplot2::ggplot(iupredDF, ggplot2::aes(x = Position))
        if (!is.null(redoxSenstitiveDF[1, 1])) {
            gg <- gg + ggplot2::geom_rect(inherit.aes = FALSE,
                    data = redoxSenstitiveDF, alpha = 0.5, fill = "#5DD39E",
                    ggplot2::aes_(xmin = ~ V1, xmax = ~ V2, ymin = 0, ymax = 1))
        }
        gg <- gg + ggplot2::geom_hline(yintercept = 0.5, size = 1, alpha = 0.5,
                                    linetype = "dotdash", color = "gray13") +
                ggplot2::geom_line(linetype = "solid",
                        ggplot2::aes(y = iupredMinus, color = "iupredMin")) +
                ggplot2::geom_line(linetype = "solid",
                        ggplot2::aes(y = iupredPlus, color = "iupredPlus")) +
                ggplot2::scale_color_manual(labels = c("Plus", "Minus"),
                        name = "Redox-Sensitive\nDisorder Prediction",
                values = c("iupredPlus" = "#BF3EFF", "iupredMin" = "#348AA7")) +
                ggplot2::labs(title = plotTitle, subtitle = plotSubtitle,
                        x = "Residue", y = "Score") + ggplot2::theme_minimal() +
                ggplot2::geom_hline(yintercept = c(0, 1), color = "gray2")
        return(gg)
    } else {
        iupredDF$redoxSensitive <- rep(FALSE, nrow(iupredDF))
        if (!is.null(redoxSenstitiveDF[1, 1])) { #overwrites if present
            senstitiveRegions <- unlist(Map(":", redoxSenstitiveDF$V1,
                                                redoxSenstitiveDF$V2))
            iupredDF$redoxSensitive <-
                seq_len(seqLength) %in% unlist(senstitiveRegions)
        }
        return(iupredDF)
    }
    }
