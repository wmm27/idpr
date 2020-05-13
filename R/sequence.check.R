#' Sequence Check Function
#'
#' This is used validate a sequence of amino acids.
#'   It can additionally be used to load an amino acid sequence.
#'   It can also be used to coerce a sequence into a specific format.
#' @param sequence amino acid sequence as a single character string
#'   or vector of single characters.
#'   It also supports a single character string that specifies
#'   the locaion of a .fasta or .fa file.
#' @param sequenceName Optional setting.
#'   \code{sequenceName = NA} for default naming.
#' @param method Required Setting.
#'   \code{method = c("Stop", "Warn", "Print")} "Stop" by default.
#'   "Stop" Reports invalid residues as an error and
#'  prevents the function from continuing.
#'   "Warn" Reports invalid residues through a warning function.
#'   "Print" Reports invalid residues as a message,
#'   but doesn't halt the function.
#'   Any invalid sequences will be reported as intended.
#'
#' @param outputType Required Setting. "String" By default.
#'   \code{outputType = c("String", "Vector", "List", "None")}
#'   Returns the sequence as a single string of amino acids.
#'   "Vector" returns the sequence as a vector of individual characters
#'   "List" returns the sequences as a list of individual characters
#'   "None" prevents the function from returning a sequence.
#' @param outputName Optional Setting. "Sequence" by default.
#' Allows for a custom value name during assignment.
#' @param nonstandardResidues Optional setting.
#'   Expands the amino acid alphabet.
#'   NA or Character vector required.
#'   Default values are "ACDEFGHIKLMNPQRSTVWY". Additional letters added here.
#'   \code{nonstandardResidues = c("O,U")}
#'   to allow Pyrrolysine (O) and Selenocysteine (U).
#' @param supressAAWarning If using nonstandardResidues,
#'    a warning will be issued.
#'   set \code{nonstandardResidues = T}
#'   to confirm addition of non-standard residues.
#' @param supressOutputMessage Set \code{supressOutputMessage = T}
#'   to prevent sequence validity message

#'
#'
#' @return A message and sequence are returned.
#'   If \code{supressOutputMessage = T}, the message is not returned.
#'   If \code{outputType = "None")}, the sequence is not returned.
#'   Otherwise, outputType will determine the format of the returned sequence.
#'   If the sequence contains an error, it will be reported
#'   based on the value of method.
#'   The Sequence will be assigned to the value "Sequence" if sequenceName
#'   is not specified. Otherwise the sequence is assigned to the value of
#'   sequenceName. This allows the sequences to be called by the user.
#' @export



sequenceCheck <- function(
  sequence,
  sequenceName = NA, #Can be character or NA (no name)

  method = "Stop", #Supports c("Stop", "Warn", "Print")
  #How the function deals with errors

  outputType = "String", #Supports c("String", "Vector", "List", "None")
  #if the sequence is to be returned as a characer string,
  #a vector of individual amino acids, or return no values

  outputName = NA,
  #Name to be the variable name that the sequence is assigned to.
  #Set as NA if outputType = "none" or if you wish to be named the default.

  nonstandardResidues = NA,
  #NA means include the typical letters,
  #expand supported residues with character vector with 1 or more residues

  supressAAWarning = F,
  #Set this to True to prevent the warning message if
  #including non standard amino acid residues

  supressOutputMessage = F
  #Set this to True to prevent a message if the sequence is valid
) {
  #-----
  #This section determines if the format of input variables are correct
  if (is.logical(supressAAWarning) == F ||
      is.logical(supressOutputMessage) == F) {
    stop("Error: supressAAWarning must be logical. Se = c(T,F)")
  }

  if (is.character(sequence) == F ||
      is.character(method) == F ||
      is.character(outputType) == F) {
    stop("Error: sequence, method, and outputType must be character vectors.
     Please check variable type.")
  }

  if ((method %in% c("Stop", "Warn", "Print")) == F) {
    stop('Error: method is not equal to a valid term.
         Set method equal to "Stop", "Warn", or "Print"')
  }



  #-----
  #This section will confirm what to do with the amino acid sequence

  if (length(sequence) == 1) {

    #this is to see if the string is a .fasta / .fa file
    if (grepl("\\.fa", sequence, ignore.case = T)) {
      sequence <- seqinr::read.fasta(file = sequence,
                                     seqtype = "AA",
                                     as.string = TRUE)
      sequence <- unlist(sequence)

      if (is.na(sequenceName)) {
    #This takes the sequence name from the fasta if  seq name is not specified
        sequenceName <- names(sequence)
      }

    }

    separatedSequence <- strsplit(sequence, "")
    names(separatedSequence) <- NULL
    separatedSequence <- unlist(separatedSequence)
  } else {
    separatedSequence <- sequence
  }

  #-----
  #This setion sets the residues which are considered valid residues

  aa <- "ACDEFGHIKLMNPQRSTVWY"
  aa <- strsplit(aa, "")
  aa <- unlist(aa)

  if (is.na(nonstandardResidues) == F) {
    aa <- c(aa, nonstandardResidues)

    if (supressAAWarning == F) {
      warningMessage <-
        paste("This includes the following non standard amino acids: ",
              nonstandardResidues,
              ". If this is an error, please set nonstandardResidues = NA . ",
              "If this is not an error, please set supressAAWarning = T",
              sep = "")
      warning(warningMessage)
    }

  }
  #-----
  #This section checks if the amino acid sequence contains invalid residues

  aaError <- F
      #Used for returning messages later. Set to False unless there is an error

  if (all(separatedSequence %in% aa) == F) {

    aaError <- T #Used for returning messages later reporting an error

    invalidResidues <- separatedSequence[!(separatedSequence %in% aa)]
    invalidResidues <- unique(invalidResidues)

    warningMessage <- paste("Protein contains the following invalid residues: ",
                            invalidResidues,
                             ". ",
                             sep = "")
    #makes the message to report what invalid residues are in the sequence

    #--- below reports the error
    if (method == "Stop") {
      stop(warningMessage)
    }

    if (method == "Warn") {
      warning(warningMessage)
    }

    if (method == "Print") {
      message(warningMessage)
    }
    #--- below reports the error
  }

  #------
  #this section creates the output

  if (is.na(outputName)) {

    if (is.na(sequenceName)) {
      outputName <- "Sequence"
    }else{
      outputName <- sequenceName
    }
  }

  if (is.na(sequenceName)) {
    sequenceName <- paste("Sequence ",
                          Sys.time(),
                          sep = "")
  }

  outputName <- gsub("[[:punct:]]", ".", outputName)

  if (outputType == "String") {
    outputSequence <- paste(separatedSequence, sep = "", collapse = "")
    names(outputSequence) <- sequenceName
    assign(outputName, outputSequence, envir = .GlobalEnv)
  }


  if (outputType == "Vector") {

    outputSequence <- separatedSequence
    assign(outputName, outputSequence, envir = .GlobalEnv)
  }

  sequenceVariableType <- is(outputSequence)[1:2] #used for output message later

  if (outputType == "List") {
     outputSequence <- separatedSequence
     outputSequence <- list(outputSequence)
     names(outputSequence) <- sequenceName
     assign(outputName, outputSequence, envir = .GlobalEnv)

    sequenceVariableType <- is(outputSequence)[1]
    #replaces the sequenceVariableType since the output in the list differs
  }

  if (supressOutputMessage == F) {

    if (aaError == F) {
      if (outputType == "None") {
        validMessage <- paste("The sequence contains no invalid residues.",
                               sep = "")
      } else {
        sequenceVariableType <- paste(sequenceVariableType, collapse = " ")
        validMessage <- paste("The sequence, named ",
                               sequenceName,
                               ", contains no invalid residues. ",
                               "The sequence was assigned to the value ",
                               outputName,
                               " as a ",
                               sequenceVariableType,
                               ".",
                               sep = "")
      }
    }

    if (aaError == T) {
      if (outputType == "None") {
        validMessage <- paste("INVALID SEQUENCE! There are invalid residues.",
                               sep = "")
      } else {
        sequenceVariableType <- paste(sequenceVariableType, collapse = " ")
        validMessage <- paste("INVALID SEQUENCE! The sequence, named ",
                               sequenceName,
                               ", DOES contains invalid residues. ",
                               "The sequence was assigned to the value ",
                               outputName,
                               " as a ",
                               sequenceVariableType,
                               ".",
                               sep = "")
      }
    }

    message(validMessage)
  }

  if (!outputType == "None") {
    return(outputSequence)
  }

}
