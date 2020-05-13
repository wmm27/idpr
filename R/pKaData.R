#' Sets of pKa values for Charged Amino Acids
#'
#' A dataset contatining the various pKa accepted values for each charged amino
#'   acid residue. N- and C-terminus values are also included. See "IPC -
#'   Isoelectric Point Calculator" Kozlowski (2016) for information on
#'   variability in pKa Data sets.
#'   \url{https://doi.org/10.1186/s13062-016-0159-9}\cr
#'   Citations are also contained in the data frame for convenience.
#'   Please cite the specific pKa set source and/or Kozlowski (2016).
#'
#' @format a data frame with 10 rows and 16 variables.
#' \describe{
#'    \item{AA}{Amino acid residues as a single letter. \cr Resudues are:
#'      Cys (C), Asp (D), Glu (E), His (H), Lys (K), Arg (R), and Tyr (Y).\cr
#'      N- and C-termini as NH2 and COOH, respectivly. "citation"
#'      in the final row}
#'    \item{EMBOSS}{pKa Dataset from
#'      \url{https://doi.org/10.1016/S0168-9525(00)02024-2}}
#'    \item{DTASelect}{pKa Dataset from
#'      \url{https://doi.org/10.1021/pr015504q}}
#'    \item{Solomons}{pKa Dataset from:
#'      \url{https://doi.org/10.1186/s13062-016-0159-9}\cr
#'          pKa Dataset original source: ISBN: 978-1-118-87576-6}
#'    \item{Sillero}{pKa Dataset from
#'      \url{https://doi.org/10.1016/0003-2697(89)90136-X}}
#'    \item{Rodwell}{pKa Dataset from
#'      \url{https://doi.org/10.1016/0003-2697(82)90611-X}}
#'    \item{Lehninger}{pKa Dataset from ISBN-13: 978-1-4641-2611-6}
#'    \item{Toseland}{pKa Dataset from
#'      \url{https://doi.org/10.1093/nar/gkj035}}
#'    \item{Thurlkill}{pKa Dataset from
#'      \url{https://doi.org/10.1110/ps.051840806}}
#'    \item{Nozaki}{pKa Dataset from:
#'      \url{https://doi.org/10.1110/ps.051840806}\cr
#'          pKa Dataset original source:
#'      \url{https://doi.org/10.1016/S0076-6879(67)11088-4}}
#'    \item{Dawson}{pKa Dataset from:
#'      \url{https://doi.org/10.1186/s13062-016-0159-9}\cr
#'          pKa Dataset original source: ISBN: 9780198552994}
#'    \item{Bjellqvist}{pKa Dataset from:
#'      \url{https://doi.org/10.1186/s13062-016-0159-9}\cr
#'          pKa Dataset original source:
#'      \url{https://doi.org/10.1002/elps.1150150171} }
#'    \item{ProMoST}{pKa Dataset from
#'      \url{https://doi.org/10.1093/nar/gkh356}}
#'    \item{IPC_protein}{pKa Dataset from
#'      \url{https://doi.org/10.1186/s13062-016-0159-9}}
#'    \item{IPC_peptide}{pKa Dataset from
#'      \url{https://doi.org/10.1186/s13062-016-0159-9}}
#'    \item{Vollhardt}{pKa Dataset from ISBN-13: 978-1-4641-2027-5}
#'  }
#' @section Additional Information:
#'  Values for NH2 and COOH are averages of values provided within the
#'  Lehringer, ProMoST, and Volhardt datasets.
#'  Lehringer and VolhardtBoth are the Seventh edition. Lehringer varies
#'  from data presented in IPC. When values could not be sourced to the original
#'  source, values were taken from Kozlowski (2016),
#'  \url{https://doi.org/10.1186/s13062-016-0159-9}. Both Kozlowski (2016) and
#'  the orignal source DOI (where avalible) or ISBN are provided within the
#'  Format section of this documention.
#'
#' @section Full Citations:
#'  Dawson, Elliott, Elliott, & Jones, 2002; Halligan et al., 2004; Kozlowski,
#'  2016; Nelson & Cox, 2017; Nozaki & Tanford, 1967; Rice, Longden, & Bleasby,
#'  2000; Rodwell, 1982; Sillero & Ribeiro, 1989; Tabb, McDonald, & Yates, 2002;
#'  TG, 1992; Thurlkill, Grimsley, Scholtz, & Pace, 2006; Toseland, McSparron,
#'  Davies, & Flower, 2006; Vollhardt & Schore, 2014)\cr
#'  Dawson, R. M. C., Elliott, D. C., Elliott, W. H., & Jones, K. M. (2002).
#'  Data for biochemical research (Vol. 3): Clarendon Press.\cr
#'  Halligan, B. D., Ruotti, V., Jin, W., Laffoon, S., Twigger, S. N., & Dratz,
#'  E. A. (2004). ProMoST (Protein Modification Screening Tool): a web-based
#'  tool for mapping protein modifications on two-dimensional gels. Nucleic
#'  Acids Research, 32(Web Server issue), W638-W644. doi:10.1093/nar/gkh356\cr
#'  Kozlowski, L. P. (2016). IPC – Isoelectric Point Calculator. Biology Direct,
#'   11(1), 55. doi:10.1186/s13062-016-0159-9\cr
#'  Nelson, D. L., & Cox, M. M. (2017). Lehninger Principles of Biochemistry
#'  (Seventh ed.). New York, NY: W. H. Freeman and Company.\cr
#'  Nozaki, Y., & Tanford, C. (1967). [84] Examination of titration behavior.
#'  In Methods in Enzymology (Vol. 11, pp. 715-734): Academic Press.\cr
#'  Rice, P., Longden, I., & Bleasby, A. (2000). EMBOSS: The European Molecular
#'  Biology Open Software Suite. Trends in Genetics, 16(6), 276-277.
#'  doi:10.1016/S0168-9525(00)02024-2\cr
#'  Rodwell, J. D. (1982). Heterogeneity of component bands in isoelectric
#'  focusing patterns. Analytical Biochemistry, 119(2), 440-449.
#'  doi:https://doi.org/10.1016/0003-2697(82)90611-X\cr
#'  Sillero, A., & Ribeiro, J. M. (1989). Isoelectric points of proteins:
#'  Theoretical determination. Analytical Biochemistry, 179(2), 319-325.
#'  doi:https://doi.org/10.1016/0003-2697(89)90136-X \cr
#'  Tabb, D. L., McDonald, W. H., & Yates, J. R. (2002). DTASelect and Contrast:
#'  Tools for Assembling and Comparing Protein Identifications from Shotgun
#'  Proteomics. Journal of Proteome Research, 1(1), 21-26. doi:10.1021/pr015504q
#'  \cr
#'  TG, S. (1992). Organic chemistry. USA: John Wiley & Sons.\cr
#'  Thurlkill, R. L., Grimsley, G. R., Scholtz, J. M., & Pace, C. N. (2006).
#'  pK values of the ionizable groups of proteins. Protein science :
#'  a publication of the Protein Society, 15(5), 1214-1218.
#'  doi:10.1110/ps.051840806\cr
#'  Toseland, C. P., McSparron, H., Davies, M. N., & Flower, D. R. (2006). PPD
#'   v1.0—an integrated, web-accessible database of experimentally determined
#'   protein pKa values. Nucleic Acids Research, 34(suppl_1), D199-D203.
#'   doi:10.1093/nar/gkj035\cr
#'  Vollhardt, P., & Schore, N. (2014). Organic Chemistry: Structure and
#'  Function (Seventh ed.). New York, NY: W. H. Freeman and Company.\cr
#' @source Kozlowski, L. P. (2016). IPC – Isoelectric Point Calculator.
#' Biology Direct, 11(1), 55.
#' \url{doi:10.1186/s13062-016-0159-9}
#' @seealso \code{\link{hendersonHasselbalch}}
"pKaData"
