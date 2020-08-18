
<!-- README.md is generated from README.Rmd. Please edit that file -->

# idpr

Overall, ‘idpr’ aims to integrate tools for the computational analysis
of intrinsically disordered proteins within R. This package is used to
identify known characteristics of IDPs within a sequence of interest
with easily reported and dynamic results. Additionally, this package
also includes tools for IDP-based sequence analysis to be used in
conjunction with other R packages.

**Please Refer to idpr-vignette.Rmd file for a detailed introduction to
the** **idpr package.**

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools") #if not already installed
devtools::install_github("wmm27/idpr")
```

## Example

This is a basic example to quickly profile your protein of interest:

``` r
library(idpr)

P53_HUMAN <- idpr:::TP53Sequences[2] #Getting a preloaded sequence from idpr
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
