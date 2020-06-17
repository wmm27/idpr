#' A Substitution Matrix for Aligning Intrinsically Disordered Proteins
#'
#' \strong{This matrix was developed and described in
#'   \href{https://doi.org/10.1142/9789812799623_0055}{Radivojac et al. (2002)}.
#'    }\cr The name "DUNMat" is taken from
#'   \href{https://doi.org/10.1038/s41598-019-52532-8}{Trivedi and
#'    Nagarajaram (2019)}. This is to keep naming consistent and distinct from
#'    other matrices named "Disorder".\cr
#'    In short: This is a substitution scoring matrix used to align proteins or
#'    regions which experience intrinsic disorder. The scores for this matrix
#'    are derived from proteins that have long regions of disorder (LDR),
#'    defined in this paper as an intrinsically disordered region (IDR) of at
#'    least 40 sequential residues. 55 protein families with LDRs were used to
#'    generate the data. Direct comparisons were not made against BLOSUM or PAM
#'    matrices within the source paper due to differences in scaling, however,
#'    when ranking its performance, it preformed the best in aligning proteins
#'    with less than 50\% sequence identity. Please see the source material,
#'    specifically, table 2, for additional information. \cr \cr
#'   \href{https://doi.org/10.1038/s41598-019-52532-8}{Trivedi and
#'    Nagarajaram (2019)} compared EDSSMat62 and DUNMat and show that DUNMat, on
#'    average, attained smaller E-values in the dataset of IDPs enriched in
#'    ordered regions, while EDSSMat62 attained smaller E-values in sets of
#'    highly disordered IDPs. Please
#'    see the referenced paper, specifically Supplementary Figure S21, for
#'    additional information and original comparison.  \cr \cr
#'    Additionally, please cite the source article when using the "DUNMat"
#'    Matrix.
#' @format A symmetrical matrix. 20x20 representing the 20 standard amino acids
#' @section Optimal Gap Parameters:
#'   These values were described in the source article and reported in Table 2.
#'   After the optimal parameters were determined, the authors further refined
#'   the gap costs.
#'   Therefore, it is recommended to use these parameters for any alignment
#'   utilizing this matrix. These were:
#'   \tabular{ccc}{
#'   DUNMat \tab Gap Open \tab Gap Extension \cr
#'   Original Optimization \tab -3 \tab -0.5 \cr
#'   Further Refinement \tab -3.2 \tab -0.1} \cr
#'   It should also be noted that a more recent work,
#'   \href{https://doi.org/10.1038/s41598-019-52532-8}{Trivedi and
#'    Nagarajaram (2019)}, determined optimal parameters based on the disordered
#'    content of query sequences, as reported in the paper's Supplementary Table
#'    S5.
#'   \tabular{ccccccc}{
#'     Matrix Name \tab Gap Open (LD) \tab Gap Extension (LD) \tab Gap Open (MD)
#'      \tab Gap Extension (MD) \tab Gap Open (HD) \tab Gap Extension (HD) \cr
#'     DUNMat \tab -6 \tab  -1 \tab  -6 \tab  -1 \tab  -16 \tab  -2}
#'    Please see the referenced paper for additional information and original
#'    reporting. Additionally, please see \code{\link{EDSSMat}}.
#'
#' @source Radivojac, P., Obradovic, Z., Brown, C. J., & Dunker, A. K. (2001).
#'   Improving sequence alignments for intrinsically disordered proteins. In
#'   Biocomputing 2002 (pp. 589-600): World Scientific.
#'   \url{https://doi.org/10.1142/9789812799623_0055}
#' @section Additional Reference:
#'   Trivedi, R., Nagarajaram, H.A. Amino acid substitution scoring
#'   matrices specific to intrinsically disordered regions in proteins. Sci
#'   Rep 9, 16380 (2019).
#'   \url{https://doi.org/10.1038/s41598-019-52532-8}
#' @family IDP-based Substitution Matrices
#' @seealso EDSSMat62

#'
"DUNMat"
