library(testthat)
library(idpr)

context("Sequence Check and Conversion")

test_that("sequenceCheck checks for accepted amino acids", {

    #--- Create vectors of amino acids

    aaStandard <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
    aaNonstandard <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                       "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

    #--- Test Accepted
    expect_message(sequenceCheck(sequence = aaStandard),
                   "The sequence contains no invalid residues.")

    #--- Test Rejected
    expect_error(sequenceCheck(sequence = aaNonstandard),
                 "Protein contains the following invalid residues: B. ")

    #--- Test changing AA
    expect_message(sequenceCheck(sequence = aaStandard,
                               nonstandardResidues = c("B"),
                               suppressAAWarning = TRUE),
                   "The sequence contains no invalid residues.")

    expect_message(sequenceCheck(sequence = aaNonstandard,
                                 nonstandardResidues = c("B"),
                                 suppressAAWarning = TRUE),
                   "The sequence contains no invalid residues.")

})


test_that("sequenceCheck accepts and returns various formats", {

    #--- Create Output
    vect_2_vect <- sequenceCheck(sequence = idpr:::HIM5Vector,
                                 outputType = "vector",
                                 suppressOutputMessage = TRUE)

    str_2_str <- sequenceCheck(sequence = idpr:::HIM5String,
                                 outputType = "string",
                                 suppressOutputMessage = TRUE)

    vect_2_str <- sequenceCheck(sequence = idpr:::HIM5Vector,
                               outputType = "string",
                               suppressOutputMessage = TRUE)

    str_2_vect <- sequenceCheck(sequence = idpr:::HIM5String,
                                 outputType = "vector",
                                 suppressOutputMessage = TRUE)
    #---

    #--- sequenceCheck removes vector name
    HIM5String_NoName <- idpr:::HIM5String
    names(HIM5String_NoName) <- NULL
    #---

    expect_equal(vect_2_vect, idpr:::HIM5Vector)
    expect_equal(str_2_str, HIM5String_NoName)
    expect_equal(str_2_vect, idpr:::HIM5Vector)
    expect_equal(vect_2_str, HIM5String_NoName)
})
