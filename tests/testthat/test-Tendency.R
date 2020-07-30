library(testthat)
library(idpr)

context("Structural Tendency Tests")

test_that("structuralTendency defaults are correctly identified", {
    P_DF <- structuralTendency("P")
    T_DF <- structuralTendency("T")
    Y_DF <- structuralTendency("Y")

    expect_equal(P_DF$Tendency,
                 "Disorder Promoting")

    expect_equal(T_DF$Tendency,
                 "Disorder Neutral")

    expect_equal(Y_DF$Tendency,
                 "Order Promoting")

})

test_that("structuralTendency matches are correctly identified after changes", {
    P_DF <- structuralTendency(sequence = "P",
                               disorderPromoting = c("T"),
                               disorderNeutral = c("Y"),
                               orderPromoting = c("P"))
    T_DF <- structuralTendency(sequence = "T",
                               disorderPromoting = c("T"),
                               disorderNeutral = c("Y"),
                               orderPromoting = c("P"))
    Y_DF <- structuralTendency(sequence = "Y",
                               disorderPromoting = c("T"),
                               disorderNeutral = c("Y"),
                               orderPromoting = c("P"))

    expect_equal(P_DF$Tendency,
                 "Order Promoting")

    expect_equal(T_DF$Tendency,
                 "Disorder Promoting")

    expect_equal(Y_DF$Tendency,
                 "Disorder Neutral")

})

test_that("structuralTendencyPlot returns expected values", {
  testDF <- structuralTendencyPlot("ACDEFGHIKLMNPQRSTVWY",
                                   graphType = "none")

  nDisorderNeutral <- table(testDF$Tendency)[1]
  names(nDisorderNeutral) <- NULL
  nDisorderPromoting <- table(testDF$Tendency)[2]
  names(nDisorderPromoting) <- NULL
  nOrderPromoting <- table(testDF$Tendency)[3]
  names(nOrderPromoting) <- NULL

  expect_equal(nDisorderNeutral,
               3)
  expect_equal(nDisorderPromoting,
               7)
  expect_equal(nOrderPromoting,
               10)
})
