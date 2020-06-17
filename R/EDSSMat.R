#' EDSSMat Disorder-based Substitution Matrices.
#'
#' The EDSSMat series of matrices were deveveloped and described in
#'   \href{https://doi.org/10.1038/s41598-019-52532-8}{Trivedi and
#'    Nagarajaram (2019)}. \cr
#'    In short: These are substitution scoring matrices
#'    used to align proteins or regions which experience intrinsic disorder.
#'    Alignment blocks, used to compute the matrix values, were composed of
#'    predicted intrinsically disordered regions. When compared to other, more
#'    frequently used substitution matrices (like BLOSUM and PAM), EDSSMat
#'    had significantly smaller E-values when aligning regions of disorder.
#'    Additionally, EDSSMat62 was shown to identify both close and distant
#'    homologs of a specific IDP while other matrices could only identify some
#'    close homologs. See the source article for additional information
#'    and for comparisons to other matrices. \cr \cr
#'    Additionally, please cite the source article when using any
#'    EDSSMat matrix.
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
#' @section Matrices:
#'   There are 7 reported EDSSMat matrices. Each vary depending on the percent
#'   identity threshold used to cluster protein sequences.
#'   EDSSMat50 clustered proteins with 50\% identity or higher,
#'   EDSSMat62 clustered proteins with 62\% identity or higher, etc. \cr
#'   \strong{See Usage Section for available matrices}
#'
#' @section Optimal Gap Parameters:
#'   These values were described in the source article, and reported in
#'   Supplemental Table S5. Therefore, it is recommended to use these parameters
#'   for any alignment utilizing the respective EDSS matrix. These were
#'   determined for 3 categories: Proteins containing Less Disorder (LD),
#'   defined as [0-20\%] disorder, Moderate Disorder (MD), defined as (20-40\%]
#'   disorder, and High Disorder (HD), defined as (40-100\%] disorder. \cr
#'   Please see the source article for additional information.
#'   \tabular{ccccccc}{
#'     Matrix Name \tab Gap Open (LD) \tab Gap Extension (LD) \tab Gap Open (MD)
#'      \tab Gap Extension (MD) \tab Gap Open (HD) \tab Gap Extension (HD) \cr
#'     EDSSMat60 \tab -7 \tab -1 \tab -6 \tab -2 \tab -14 \tab -3 \cr
#'     EDSSMat62 \tab -8 \tab -1 \tab -5 \tab -2 \tab -19 \tab -2 \cr
#'     EDSSMat70 \tab -7 \tab -1 \tab -5 \tab -2 \tab -19 \tab -2 \cr
#'     EDSSMat75 \tab -8 \tab -1 \tab -5 \tab -2 \tab -19 \tab -2 \cr
#'     EDSSMat80 \tab -7 \tab -1 \tab -5 \tab -2 \tab -15 \tab -3 \cr
#'     EDSSMat90 \tab -7 \tab -1 \tab -5 \tab -2 \tab -19 \tab -2 }
#'
#' @source Trivedi, R., Nagarajaram, H.A. Amino acid substitution scoring
#'   matrices specific to intrinsically disordered regions in proteins. Sci
#'   Rep 9, 16380 (2019).
#'   \url{https://doi.org/10.1038/s41598-019-52532-8}
#' @family IDP-based Substitution Matrices
#' @name EDSSMat

NULL


#' @rdname EDSSMat
"EDSSMat50"

#' @rdname EDSSMat
"EDSSMat60"

#' @rdname EDSSMat
"EDSSMat62"

#' @rdname EDSSMat
"EDSSMat70"

#' @rdname EDSSMat
"EDSSMat75"

#' @rdname EDSSMat
"EDSSMat80"

#' @rdname EDSSMat
"EDSSMat90"
