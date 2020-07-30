library(testthat)
library(idpr)

context("Charge Tests")

test_that("netCharge returns the same values for various input types", {
    expect_equal(netCharge(idpr:::HIM5String),
                 netCharge(idpr:::HIM5Vector))

    expect_equal(netCharge(idpr:::HIM5String,
                           pH = 5),
                 netCharge(idpr:::HIM5Vector,
                           pH = 5))

    expect_equal(netCharge(idpr:::HIM5String,
                           pKaSet = "EMBOSS"),
                 netCharge(idpr:::HIM5Vector,
                           pKaSet = "EMBOSS"))
})

test_that("chargeCalculationGlobal termini calculates correctly", {
    withTerminiDF <- chargeCalculationGlobal(idpr:::HIM5String,
                                             plotResults = FALSE,
                                             includeTermini = TRUE)
    withoutTerminiDF <- chargeCalculationGlobal(idpr:::HIM5String,
                                                plotResults = FALSE,
                                                includeTermini = FALSE)
    expect_equal(withTerminiDF$Charge[1],
                 hendersonHasselbalch(pKa = 9.0939999999999994,
                                      pH = 7,
                                      residue = "NH2"))

    expect_equal(withoutTerminiDF$Charge[1],
                 0)

    expect_equal(withTerminiDF$Charge[3],
                 hendersonHasselbalch(pKa = 11.84,
                                      pH = 7,
                                      residue = "R"))

    expect_equal(withoutTerminiDF$Charge[3],
                 hendersonHasselbalch(pKa = 11.84,
                                      pH = 7,
                                      residue = "R"))
})

test_that("chargeCalculationLocal window size change", {
    window1DF <- chargeCalculationLocal(idpr:::HIM5String,
                                        plotResults = FALSE,
                                        window = 1)

    window3DF <- chargeCalculationLocal(idpr:::HIM5String,
                                        plotResults = FALSE,
                                        window = 3)

    window7DF <- chargeCalculationLocal(idpr:::HIM5String,
                                        plotResults = FALSE,
                                        window = 7)

    window11DF <- chargeCalculationLocal(idpr:::HIM5String,
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
