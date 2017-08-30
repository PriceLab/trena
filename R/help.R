#' Inferring Transcriptional Regulation with TReNA
#'
#' `trena` provides a framework for using gene expression data to infer relationships
#' between a target gene and a set of transcription factors. It does so using a
#' several classes and their associated methods, briefly documented below
#'
#' #'
#' Solver Class Objects 
#'
#'  The \code{\link{Solver}} class is a base class used within `trena`. A particular
#' \code{\link{Solver}} object also contains the name of the selected solver and
#' dispatches the correct feature selection method when \code{run} is called on the
#' object. It is inherited by all the following subclasses, representing the
#' different feature selection methods: \code{\link{BayesSpikeSolver}},
#' \code{\link{EnsembleSolver}}, \code{\link{LassoPVSolver}},
#' \code{\link{LassoSolver}}, \code{\link{PearsonSolver}},
#' \code{\link{RandomForestSolver}}, \code{\link{RidgeSolver}},
#' \code{\link{SpearmanSolver}}, \code{\link{SqrtLassoSolver}}.
#' 
#' CandidateFilter Class Objects
#'
#' The \code{\link{CandidateFilter}} class is a base class that is generally used to filter
#' the transcription factors in the expression matrix to obtain a set of candidate
#' regulators. Filtering method depends on the filter type chosen; there are currently the
#' following subclasses: \code{\link{FootprintFilter}}, \code{\link{HumanDHSFilter}},
#' \code{\link{GeneOntologyFilter}}, and \code{\link{VarianceFilter}}. The filters are
#' applied using the \code{\link{getCandidates}} method on a given
#' \code{\link{CandidateFilter}} object. 
#'
#' FootprintFinder Class Objects
#'
#' The  \code{\link{FootprintFinder}} class is designed to allow extraction
#' of gene footprinting information from existing PostgreSQL or SQLite
#' databases. In standard use of the `trena` package, it is used solely by
#' the \code{\link{getCandidates}} method for a \code{\link{FootprintFilter}}
#' object. However, a  \code{\link{FootprintFinder}} object has many more
#' available methods that allow it to extract information more flexibly.
#' 
#' @name trena-package
#' @aliases trena-package
#'
#' @docType package

NULL
