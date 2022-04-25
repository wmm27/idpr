
<!-- README.md is generated from README.Rmd. Please edit that file -->

# idpr

Overall, ‘idpr’ aims to integrate tools for the computational analysis
of intrinsically disordered proteins within R. This package is used to
identify known characteristics of IDPs within a sequence of interest
with easily reported and dynamic results. Additionally, this package
also includes tools for IDP-based sequence analysis to be used in
conjunction with other R packages. See our recently published
peer-reviewed publication in [PLOS ONE
(https://doi.org/10.1371/journal.pone.0266929)](https://doi.org/10.1371/journal.pone.0266929)

**Please Refer to idpr-vignette.Rmd file for a detailed introduction to
the** **idpr package.**

Links to the vignettes found at the [Bioconductor landing page
(here)](https://doi.org/doi:10.18129/B9.bioc.idpr) or

## Installation

You can install the stable release version version from
[Bioconductor](https://doi.org/doi:10.18129/B9.bioc.idpr) with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("idpr")
```

Additionally, you can install the development version from
[Bioconductor](https://bioconductor.org/packages/devel/bioc/html/idpr.html)
with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')
```

Or you can install the most recent development version from
[GitHub](https://github.com/wmm27/idpr) with:

``` r
# install.packages("devtools") #if not already installed
devtools::install_github("wmm27/idpr")
```

## Example

This is a basic example to quickly profile your protein of interest:

``` r
library(idpr)

P53_HUMAN <- TP53Sequences[2] #Getting a pre-loaded sequence from idpr
print(P53_HUMAN)
#>                                                                                                                                                                                                                                                                                                                                                                                            P04637|P53_HUMAN 
#> "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"

P53_ID <- "P04637" #Human TP53 UniProt ID

#Generates the IDP Profile:
idprofile(sequence = P53_HUMAN, 
          uniprotAccession = P53_ID, 
          proteinName = "TP53 Human", #Optional Argument
          window = 11, #Optional Argument
          pKaSet = "Lehninger", #Optional Argument
          iupredType = "redox" #Optional Argument
)
#> [[1]]
```

<img src="man/figures/README-example-1.png" width="75%" />

    #> 
    #> [[2]]

<img src="man/figures/README-example-2.png" width="75%" />

    #> 
    #> [[3]]

<img src="man/figures/README-example-3.png" width="75%" />

    #> 
    #> [[4]]

<img src="man/figures/README-example-4.png" width="75%" />

    #> 
    #> [[5]]

<img src="man/figures/README-example-5.png" width="75%" />

    #> 
    #> [[6]]

<img src="man/figures/README-example-6.png" width="75%" />

**Please Refer to idpr-vignette.Rmd file for a detailed introduction to
the** **idpr package.** [Link to the Vignette
(here)](https://bioconductor.org/packages/release/bioc/vignettes/idpr/inst/doc/idpr-vignette.html)

## Appendix

For use and details on ‘idpr’, see our peer-reviewed article published
in [PLOS ONE
(https://doi.org/10.1371/journal.pone.0266929)](https://doi.org/10.1371/journal.pone.0266929)

### Package citation

``` r
citation("idpr")
#> 
#> To cite idpr in publications, use the citation below and other
#> function-specific sources found in the idpr package documentation:
#> 
#>   McFadden, W. M., and Yanowitz, J. L. (2022). idpr: A package for
#>   profiling and analyzing Intrinsically Disordered Proteins in R. PLOS
#>   ONE, 17(4), e0266929. doi:10.1371/journal.pone.0266929
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {idpr: A package for profiling and analyzing Intrinsically Disordered Proteins in R},
#>     author = {William M McFadden and Judith L Yanowitz},
#>     journal = {PLOS ONE},
#>     publisher = {Public Library of Science},
#>     year = {2022},
#>     volume = {17},
#>     number = {4},
#>     pages = {e0266929},
#>     doi = {10.1371/journal.pone.0266929},
#>     url = {https://bioconductor.org/packages/idpr/},
#>   }
```

### Function citations

-   Bálint Mészáros, Gábor Erdős, Zsuzsanna Dosztányi, IUPred2A:
    context-dependent prediction of protein disorder as a function of
    redox state and protein binding, Nucleic Acids Research, Volume 46,
    Issue W1, 2 July 2018, Pages W329–W337,
    <https://doi.org/10.1093/nar/gky384>
-   Erdős, G., & Dosztányi, Z. (2020). Analyzing protein disorder with
    IUPred2A. Current Protocols in Bioinformatics, 70, e99.
    <https://doi.org/10.1002/cpbi.99>
-   Kozlowski, L. P. (2016). IPC – Isoelectric Point Calculator. Biology
    Direct, 11(1), 55. <https://doi.org/10.1186/s13062-016-0159-9>
-   Kyte, J., & Doolittle, R. F. (1982). A simple method for displaying
    the hydropathic character of a protein. Journal of molecular
    biology, 157(1), 105-132.
-   Nelson, D. L., & Cox, M. M. (2017). Lehninger Principles of
    Biochemistry (Seventh ed.). New York, NY: W. H. Freeman and Company.
-   Prilusky, J., Felder, C. E., et al. (2005). FoldIndex: a simple tool
    to predict whether a given protein sequence is intrinsically
    unfolded. Bioinformatics, 21(16), 3435-3438.
-   Uversky, V. N. (2016). Paradoxes and wonders of intrinsic disorder:
    Complexity of simplicity. Intrinsically Disordered Proteins, 4(1),
    e1135015. <https://doi.org/10.1080/21690707.2015.1135015>
-   Uversky, V. N. (2013). A decade and a half of protein intrinsic
    disorder: Biology still waits for physics. Protein Science, 22(6),
    693-724. <doi:10.1002/pro.2261>
-   Uversky, V. N., Gillespie, J. R., & Fink, A. L. (2000). Why are
    “natively unfolded” proteins unstructured under physiologic
    conditions?. Proteins: structure, function, and bioinformatics,
    41(3), 415-427.
    <https://doi.org/10.1002/1097-0134(20001115)41:3>\<415::AID-PROT130\>3.0.CO;2-7

### Additional Information

``` r
Sys.time()
#> [1] "2022-04-25 12:20:36 EDT"
Sys.Date()
#> [1] "2022-04-25"
R.version
#>                _                           
#> platform       x86_64-apple-darwin17.0     
#> arch           x86_64                      
#> os             darwin17.0                  
#> system         x86_64, darwin17.0          
#> status                                     
#> major          4                           
#> minor          1.3                         
#> year           2022                        
#> month          03                          
#> day            10                          
#> svn rev        81868                       
#> language       R                           
#> version.string R version 4.1.3 (2022-03-10)
#> nickname       One Push-Up
```
