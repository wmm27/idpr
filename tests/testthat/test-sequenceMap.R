library(testthat)
library(idpr)

context("sequenceMapCoordinates Tests")


test_that("sequenceMapCoordinates nbResdues argument functions properly", {

    standardDF <- sequenceMapCoordinates("ACDEFGHIKLMNPQRSTVWY")

    #test if all in standardDF row are 1, due to nbResidues > seqLength
    expect_true(all(standardDF$row == 1))
    expect_false(all(standardDF$col == 1))


    residue1DF <- sequenceMapCoordinates("ACDEFGHIKLMNPQRSTVWY",
                                         nbResidues = 1)
    #test if all in residue1DF col are 1, due to nbResidues = 1
    expect_true(all(residue1DF$col == 1))
    expect_false(all(residue1DF$row == 1))

})
