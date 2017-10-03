#----------------------------------------------------------------------------------------------------
#' @name Solver-class
#' @rdname Solver-class
#' @aliases Solver
#' 
#' @import methods

.Solver <- setClass ("Solver",
                     slots = c(mtx.assay="matrix",
                               targetGene="character",
                               candidateRegulators="character",
                               quiet="logical",
                               state="environment")
                     )

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
#' Retrieve the assay matrix of gene expression data from a Solver object
#' 
#' @rdname getAssayData
#' @aliases getAssayData
#' 
#' @param obj An object of class Solver
#'
#' @export
#' 
#' @return The assay matrix of gene expression data associated with a Solver object
#'
#' @examples
#' # Create a Solver object using the included Alzheimer's data and retrieve the matrix
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' targetGene <- "MEF2C"
#' candidateRegulators <- setdiff(rownames(mtx.sub), targetGene)
#' solver <- Solver(mtx.sub, targetGene, candidateRegulators)
#' mtx <- getAssayData(solver)
#' 
setGeneric("getAssayData",    signature="obj", function(obj) standardGeneric ("getAssayData"))

#' @export
setGeneric("run",             signature="obj", function(obj) standardGeneric ("run"))

#' @export
setGeneric("rescalePredictorWeights",
           signature="obj", function(obj, rawValue.min, rawValue.max, rawValues) standardGeneric ("rescalePredictorWeights"))

#' Retrieve the target gene from a Solver object
#'
#' @rdname getTarget
#' @aliases getTarget
#'
#' @param obj An object of class Solver
#'
#' @return The target gene associated with a Solver object
#'
#' @examples
#' # Create a Solver object using the included Alzheimer's data and retrieve the target gene
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' targetGene <- "MEF2C"
#' candidateRegulators <- setdiff(rownames(mtx.sub), targetGene)
#' solver <- Solver(mtx.sub, targetGene, candidateRegulators)
#' target <- getTarget(solver) 

#' @export
setGeneric("getTarget", signature = "obj", function(obj) standardGeneric("getTarget"))

#' Retrieve the candiate regulators from a Solver object
#'
#' @rdname getRegulators
#' @aliases getRegulators
#'
#' @param obj An object of class Solver
#'
#' @return The candidate regulators associated with a Solver object
#'
#' @examples
#' # Create a Solver object using the included Alzheimer's data and retrieve the regulators
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' targetGene <- "MEF2C"
#' candidateRegulators <- setdiff(rownames(mtx.sub), targetGene)
#' solver <- Solver(mtx.sub, targetGene, candidateRegulators)
#' regs <- getRegulators(solver) 

#' @export
setGeneric("getRegulators", signature = "obj", function(obj) standardGeneric("getRegulators"))


#----------------------------------------------------------------------------------------------------
#' Define an object of class Solver
#'
#' @description
#' The Solver class is a base class that governs the different solvers available in \code{trena}.
#' It is rarely called by itself; rather, interaction with a particular solver object
#' is achieved using a specific solver type.
#' 
#' @rdname Solver-class
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param targetGene A designated target gene that should be part of the mtx.assay data
#' @param candidateRegulators The designated set of transcription factors that could be associated
#' @param quiet A logical indicating whether or not the Solver object should print output
#'
#' @export
#'
#' @return An object of the Solver class
#'
#' @examples
#'#' # Create a Solver object using the included Alzheimer's data
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' targetGene <- "MEF2C"
#' candidateRegulators <- setdiff(rownames(mtx.sub), targetGene)
#' solver <- Solver(mtx.sub, targetGene, candidateRegulators) # Create a simple Solver object with default options
#'
#' @seealso \code{\link{getAssayData}}, \code{\link{getTarget}}, \code{\link{getRegulators}}
#'
#' @family Solver class objects

Solver <- function(mtx.assay=matrix(), targetGene, candidateRegulators, quiet=TRUE)
{
    # If a matrix is supplied, check the distribution to see if it's too big
    # Also check to make sure target gene is well-enough expressed
    if(!is.na(max(mtx.assay))){
        mtx.ratio <- (max(mtx.assay) - stats::quantile(mtx.assay,0.75))/(stats::quantile(mtx.assay,0.75) - stats::median(mtx.assay))
        if(mtx.ratio > 1000){
            warning("Assay matrix may contain highly skewed data; consider transforming your matrix.")
        }
        
        if(rowMeans(mtx.assay)[targetGene] < stats::quantile(rowMeans(mtx.assay), probs = 0.1)){
            warning("Target gene mean expression is in the bottom 10% of all genes in the assay matrix")
        } 
    }

    # Check to make sure the candidate regulators and target gene aren't empty;
    # If they are, send up a warning
    if(length(targetGene) == 0) {
        warning("No target gene supplied; please supply a target gene to avoid errors")
    }
    if(length(candidateRegulators) == 0) {
        warning("No regulators supplied; please supply regulators to avoid errors")
        }
    
    env <- new.env(parent=emptyenv())
    .Solver(mtx.assay=mtx.assay,
            targetGene = targetGene,
            candidateRegulators = candidateRegulators,
            quiet=quiet,
            state=env)
    
} # Solver, the constructor
#----------------------------------------------------------------------------------------------------
setMethod("getAssayData", "Solver",
          
          function (obj){
              obj@mtx.assay
          })
#----------------------------------------------------------------------------------------------------
setMethod("getTarget", "Solver",
          
          function (obj){
              obj@targetGene
          })
#----------------------------------------------------------------------------------------------------
setMethod("getRegulators", "Solver",
          
          function (obj){
              obj@candidateRegulators
          })
#----------------------------------------------------------------------------------------------------
#' Rescale the Predictor Weights
#'
#' Solvers such as LASSO penalize predictors on a scale of 1 (full weight) to infinity (zero weight).
#' With the \code{rescalePredictorWeights} method, incoming raw values can be scaled between a possibly
#' theoretical minimum and maximum value.
#'
#' @rdname rescalePredictorWeights
#' @aliases rescalePredictorWeights
#'
#' @param obj An object of the Solver class
#' @param rawValue.min The minimum value of the raw expression values
#' @param rawValue.max The maximum value of the raw expression values
#' @param rawValues A matrix of raw expression values
#'
#' @export
#'
#' @return A matrix of the raw values re-scaled using the minimum and maximum values
#'
#' @examples
#' # Create a LassoSolver object using the included Alzheimer's data and rescale the predictors
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' targetGene <- "MEF2C"
#' candidateRegulators <- setdiff(rownames(mtx.sub), targetGene)
#' ls <- LassoSolver(mtx.sub, targetGene, candidateRegulators)
#' raw.values <- c(241, 4739, 9854, 22215, 658334)
#' cooked.values <- rescalePredictorWeights(ls, rawValue.min = 1, rawValue.max = 1000000, raw.values)

setMethod("rescalePredictorWeights", "Solver",
          
          function(obj, rawValue.min, rawValue.max, rawValues){
              1 - ((rawValues-rawValue.min)/(rawValue.max-rawValue.min))
          })
#----------------------------------------------------------------------------------------------------
