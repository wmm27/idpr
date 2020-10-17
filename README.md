
<!-- README.md is generated from README.Rmd. Please edit that file -->

# idpr

Overall, ‘idpr’ aims to integrate tools for the computational analysis
of intrinsically disordered proteins within R. This package is used to
identify known characteristics of IDPs within a sequence of interest
with easily reported and dynamic results. Additionally, this package
also includes tools for IDP-based sequence analysis to be used in
conjunction with other R packages.

**Please Refer to idpr-vignette.Rmd file for a detailed introduction to
the** **idpr package.** Links to the vignettes found at the
[Bioconductor landing page](https://doi.org/doi:10.18129/B9.bioc.idpr)

## Installation

You can install the development version from
[Bioconductor](https://doi.org/doi:10.18129/B9.bioc.idpr) with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("idpr")
```

Or you can install the development version from
[GitHub](https://github.com/wmm27/idpr) with:

``` r
# install.packages("devtools") #if not already installed
devtools::install_github("wmm27/idpr")
```

## Example

This is a basic example to quickly profile your protein of interest:

``` r
library(idpr)

P53_HUMAN <- TP53Sequences[2] #Getting a preloaded sequence from idpr
print(P53_HUMAN)
#>                                                                                                                                                                                                                                                                                                                                                                                            P04637|P53_HUMAN 
#> "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"

P53_ID <- "P04637" #Human TP53 UniProt ID

idprofile(sequence = P53_HUMAN, #Generates the Profi
          uniprotAccession = P53_ID)
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

**Please Refer to idpr-vignette.Rmd file for a detailed introduction to
the** **idpr package.**

## Appendix

### Package citation

``` r
citation("idpr")
#> 
#> To cite package 'idpr' in publications use:
#> 
#>   William McFadden and Judith Yanowitz (2020). idpr: Profiling and
#>   Analyzing Intrinsically Disordered Proteins in R. R package version
#>   0.99.25.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {idpr: Profiling and Analyzing Intrinsically Disordered Proteins in R},
#>     author = {William McFadden and Judith Yanowitz},
#>     year = {2020},
#>     note = {R package version 0.99.25},
#>   }
```

### Additional Information

``` r
Sys.time()
#> [1] "2020-10-17 14:40:21 EDT"
Sys.Date()
#> [1] "2020-10-17"
R.version
#>                _                           
#> platform       x86_64-apple-darwin17.0     
#> arch           x86_64                      
#> os             darwin17.0                  
#> system         x86_64, darwin17.0          
#> status                                     
#> major          4                           
#> minor          0.2                         
#> year           2020                        
#> month          06                          
#> day            22                          
#> svn rev        78730                       
#> language       R                           
#> version.string R version 4.0.2 (2020-06-22)
#> nickname       Taking Off Again
```
