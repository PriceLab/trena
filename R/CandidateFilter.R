#' @import methods
#'
#' @name CandidateFilter-class
#' @rdname CandidateFilter-class
#' @aliases CandidateFilter
#'
#' @slot quiet A logical denoting whether or not the CandidateFilter object should print output
#'
#' @family Filtering objects

.CandidateFilter <- setClass("CandidateFilter",
                    slots = c(quiet = "logical"
			      )
                    )
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
#' Get candidate genes using a CandidateFilter object
#'
#' @rdname getCandidates
#' @aliases getCandidates
#'
#' @param obj An object of a CandidateFilter class
#'
#' @return A vector containing genes from the assay matrix that are selected by the filter
#'
#' @family getCandidate Methods
#' @export
setGeneric("getCandidates", signature="obj", function(obj, ...) standardGeneric("getCandidates"))

#----------------------------------------------------------------------------------------------------
#' CandidateFilter
#'
#' A CandidateFilter is an S4 class to represent a gene candidate filter. These filters can employ
#' a variety of methods to reduce the number of transcription factors used as predictors for solving
#' a Solver object.
#'
#' @rdname CandidateFilter-class
#' @aliases CandidateFilter
#'
#' @param quiet A logical denoting whether or not the CandidateFilter object should print output
#'
#' @return An object of the Candidate filter class
#'
#' @seealso \code{\link{getCandidates}}
#'
#' @export
#'
#' @examples
#' # Create an empty candidate filter
#' candidate.filter <- CandidateFilter(quiet=TRUE)

CandidateFilter <- function(quiet = TRUE)
{
    .CandidateFilter(quiet = quiet)

} # CandidateFilter, the constructor
#----------------------------------------------------------------------------------------------------
