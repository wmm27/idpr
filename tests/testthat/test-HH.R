library(testthat)
library(idpr)

context("Henderson-Hasselbalch Tests")


test_that("when pKa == pH, the [acid] = [conjugate base]", {
  expect_equal(hendersonHasselbalch(pKa = 7, pH = 7, "acid"),
               -0.5)
  expect_equal(hendersonHasselbalch(pKa = 7, pH = 7, "base"),
               0.5)
})

test_that("when pKa >>> pH, the [acid] >>> [conjugate base]", {
  expect_equal(round(hendersonHasselbalch(pKa = 12, pH = 2, "acid"), 2),
               0)
  expect_equal(round(hendersonHasselbalch(pKa = 12, pH = 2, "base"), 2),
               1)
})

test_that("when pKa >>> pH, the [acid] >>> [conjugate base]", {
  expect_equal(round(hendersonHasselbalch(pKa = 2, pH = 12, "acid"), 2),
               -1)
  expect_equal(round(hendersonHasselbalch(pKa = 2, pH = 12, "base"), 2),
               0)
})

test_that("HH termini residues are equal", {
    expect_equal(hendersonHasselbalch(pKa = 2, residue = "COOH"),
                 hendersonHasselbalch(pKa = 2, residue = "COO"))

    expect_equal(hendersonHasselbalch(pKa = 2, residue = "NH3"),
                 hendersonHasselbalch(pKa = 2, residue = "NH2"))
})
