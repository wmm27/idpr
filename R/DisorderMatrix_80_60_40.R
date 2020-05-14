#' Disorder-based Substitution Matrices.
#'
#' \strong{The Disorder40, Disorder60, and Disorder85 Matrices were
#'   deveveloped and described in
#'   \href{https://doi.org/10.1093/molbev/msp277}{Brown et al. (2009)}.} \cr
#'    In short: There are substitution scoring matrices used to align proteins
#'    or regions which experience intrinsic disorder. The matrices were
#'    calculated using pairwise sequence alignments of protein families
#'    which here identified from 287 experimentally confirmed Intrinsially
#'    Disordered Proteins (IDPs). The IDPs contained at least 30 sequential
#'    residues of intrinsic disorder and protein families were found using
#'    BLAST.\cr There was not a comprehensive comparison to other frequently
#'    used substitution matrices (like BLOSUM and PAM) in terms of improving IDP
#'    sequence alignments. The authors note that the purpose of these
#'    matrices were to compare evolutionary characteristics of disordered and
#'    ordered proteins.  Please see the source material for additional
#'    information.\cr \cr
#'   \href{https://doi.org/10.1038/s41598-019-52532-8}{Trivedi and
#'    Nagarajaram (2019)} compared EDSSMat62 against all three Disordered
#'    Matrices. Disorder40 and Disorder85 attain lower E-values for
#'    highly disordered proteins, on average, when compared to EDSSMat62.
#'    EDSSMat62 attained lower E-values when compared to Disorder60 for aligning
#'    highly disordered proteins.
#'    EDSSMat62 preforms better than all three Disorder matrices for IDPs
#'    enriched in ordered regions. Please
#'    see the referenced paper, specifically Supplementary Figures S18-20, for
#'    additional information and original comparison.  \cr \cr
#'    Additionally, please cite the source article when using Disorder40,
#'    Disorder60, or Disorder85.
#'
#' @format All matrices are symmetric. 24 residues are represented:
#'   \itemize{
#'     \item Each of the standard 20 standard amino acids
#'     \item Four ambiguous residues:
#'        \itemize{
#'          \item B: Asparagine or Aspartic Acid (Asx)
#'          \item Z: Glutamine or Glutamic Acid (Glx)
#'          \item X: Unspecified or unknown amino acid
#'          \item *: Stop
#'        }
#'     }
#'
#' @section Optimal Gap Parameters:
#'  As mentioned in the Description, the intended use of these matrices were
#'    not to improve sequence alignments. Therefore no gap penalty values are
#'    provided.
#'
#'   It should also be noted that a more recent work,
#'   \href{https://doi.org/10.1038/s41598-019-52532-8}{Trivedi and
#'    Nagarajaram (2019)}, determined optimal parameters based on the disordered
#'    content of query sequences, as reported in the paper's Supplementary Table
#'    S5.
#'   \tabular{ccccccc}{
#'     Matrix Name \tab Gap Open (LD) \tab Gap Extension (LD) \tab Gap Open (MD)
#'      \tab Gap Extension (MD) \tab Gap Open (HD) \tab Gap Extension (HD) \cr
#'  Disorder40 \tab -20 \tab -1 \tab -7 \tab -1 \tab -7 \tab -1 \cr
#'  Disorder60 \tab -20 \tab -1 \tab -16 \tab -1 \tab -11 \tab -2 \cr
#'  Disorder85 \tab -20 \tab -1 \tab -16 \tab -1 \tab -7 \tab -2}
#'    Please see the referenced paper for additional information and original
#'    reporting. Additionally, please see \code{\link{EDSSMat}}.
#' @source Brown, C. J., Johnson, A. K., & Daughdrill, G. W. (2009).
#'   Comparing Models of Evolution for Ordered and Disordered Proteins.
#'   Molecular Biology and Evolution, 27(3), 609-621.
#'   \href{https://doi.org/10.1093/molbev/msp277}{doi:10.1093/molbev/msp277}
#' @section Additional Reference:
#'   Trivedi, R., Nagarajaram, H.A. Amino acid substitution scoring
#'   matrices specific to intrinsically disordered regions in proteins. Sci
#'   Rep 9, 16380 (2019).
#'   \url{https://doi.org/10.1038/s41598-019-52532-8}
#' @family IDP-based Substitution Matrices
#' @seealso EDSSMat62
#' @name DisorderMat
NULL


#' @rdname DisorderMat
"Disorder40"

#' @rdname DisorderMat
"Disorder60"

#' @rdname DisorderMat
"Disorder85"


