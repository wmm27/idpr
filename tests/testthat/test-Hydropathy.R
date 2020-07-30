library(testthat)
library(idpr)

context("Hydropathy Tests")

test_that("scaledHydropathyLocal window size change", {
    window1DF <- scaledHydropathyLocal(idpr:::HIM5String,
                                        plotResults = FALSE,
                                        window = 1)

    window3DF <- scaledHydropathyLocal(idpr:::HIM5String,
                                        plotResults = FALSE,
                                        window = 3)

    window7DF <- scaledHydropathyLocal(idpr:::HIM5String,
                                        plotResults = FALSE,
                                        window = 7)

    window11DF <- scaledHydropathyLocal(idpr:::HIM5String,
                                         plotResults = FALSE,
                                         window = 11)

    expect_equal(nchar(window1DF$Window[1]),
                 1)
    expect_equal(nchar(window3DF$Window[1]),
                 3)
    expect_equal(nchar(window7DF$Window[1]),
                 7)
    expect_equal(nchar(window11DF$Window[1]),
                 11)
})

test_that("meanScaledHydropathy returns expected values", {
    expect_equal(meanScaledHydropathy("R"),
                 0)
    expect_equal(meanScaledHydropathy("I"),
                 1)
    expect_equal(meanScaledHydropathy("A"),
                 idpr:::KDNorm$V2[1])
})
