#' Human Cellular tumor antigen p53 and homologs
#'
#' This is a vector of sequences as character strings. This contains the
#'    amino acid sequence of Human Cellular tumor antigen p53
#'    (UniProt ID: P04637) and sequences of several homologous sequences.
#'    These sequences in TP53Sequences selected due to their highly similar
#'    identity on  UniProt (
#'    \href{https://doi.org/10.1093/nar/gky1049}{The UniProt Consortium, 2019}
#'    ).
#'
#' @section UniProt IDs:
#'   \itemize{
#'     \item P02340
#'     \item P04637
#'     \item P10361
#'     \item Q29537
#'     \item Q00366
#'     \item O09185
#'     \item Q9TTA1
#'     \item Q95330
#'     \item A0A2I2Y7Z8*
#'     }
#'  * The Gorilla p53 sequence is not within the TP53Sequences vector, and is
#'  its own object named GorillaTP53. This is because the Gorilla p53 is
#'  unreviewed on UniProt, so we chose to exclude it from the list of otherwise
#'  SwissProt sequences.
#' @family Example Sequences
#' @source UniProt Consortium. (2019).
#'  UniProt: a worldwide hub of protein knowledge.
#'  Nucleic acids research, 47(D1), D506-D515.
#'  https://doi.org/10.1093/nar/gky1049
#' @name TP53Sequences
NULL

#' @rdname TP53Sequences
"TP53Sequences"

#' @rdname TP53Sequences
"GorillaTP53"
