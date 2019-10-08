#----------------------------------------------------------------------------------------------------
#' An S4 class to represent a XGBoost solver
#'
#' @import xgboost
#' @include Solver.R
#' @import methods
#'
#' @name XGBoostSolver-class
#'
.XGBoostSolver <- setClass("XGBoostSolver", contains = "Solver")
#----------------------------------------------------------------------------------------------------
#' Create a Solver class using gradient boosting (a regression technique) and the XGBoost library
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param targetGene A designated target gene that should be part of the mtx.assay data
#' @param candidateRegulators The designated set of transcription factors that could be associated
#' with the target gene
#' @param quiet A logical denoting whether or not the solver should print output
#'
#' @return A Solver class object with XGBoost Importances (Gain) as the solver
#'
#' @seealso  \code{\link{solve.XGBoost}}, \code{\link{getAssayData}}
#'
#' @family Solver class objects
#'
#' @export
#'
#' @examples
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' XGBoost.solver <- XGBoostSolver(mtx.sub, target.gene, tfs)

XGBoostSolver <- function(mtx.assay = matrix(), targetGene, candidateRegulators, quiet=TRUE)
{
    # Remove the targetGene from candidateRegulators
    if(any(grepl(targetGene, candidateRegulators)))
        candidateRegulators <- candidateRegulators[-grep(targetGene, candidateRegulators)]

    # Check to make sure the matrix contains some of the candidates
    candidateRegulators <- intersect(candidateRegulators, rownames(mtx.assay))
    stopifnot(length(candidateRegulators) > 0)

    obj <- .XGBoostSolver(Solver(mtx.assay=mtx.assay,
                                  quiet=quiet,
                                  targetGene=targetGene,
                                  candidateRegulators=candidateRegulators))

    # Send a warning if there's a row of zeros
    if(!is.na(max(mtx.assay)) & any(rowSums(mtx.assay) == 0))
        warning("One or more gene has zero expression; this may yield 'NA' results and warnings when using XGBoost correlations")

    obj

} #XGBoostSolver, the constructor
#----------------------------------------------------------------------------------------------------
#' Show the XGBoost Solver
#'
#' @rdname show.XGBoostSolver
#' @aliases show.XGBoostSolver
#'
#' @param object An object of the class XGBoostSolver
#'
#' @return A truncated view of the supplied object
#'
#' @examples
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' XGBoost.solver <- XGBoostSolver(mtx.sub, target.gene, tfs)
#' show(XGBoost.solver)

setMethod('show', 'XGBoostSolver',

          function(object) {
              regulator.count <- length(getRegulators(object))
              if(regulator.count > 10){
                  regulatorString <- paste(getRegulators(object)[1:10], collapse=",")
                  regulatorString <- sprintf("%s...", regulatorString);
              }
              else
                  regulatorString <- paste(getRegulators(object), collapse=",")

              msg = sprintf("XGBoostSolver with mtx.assay (%d, %d), targetGene %s, %d candidate regulators %s",
                            nrow(getAssayData(object)), ncol(getAssayData(object)),
                            getTarget(object), regulator.count, regulatorString)
              cat (msg, '\n', sep='')
          })
#----------------------------------------------------------------------------------------------------
#' Run the XGBoost Solver
#'
#' @rdname solve.XGBoost
#' @aliases run.XGBoostSolver solve.XGBoost
#'
#' @description Given a TReNA object with XGBoost as the solver, use the \code{\link{cor}}
#' function with \code{method = "XGBoost"} to esimate importances for each transcription factor
#' as a predictor of the target gene's expression level.
#'
#' @param obj An object of class XGBoostSolver
#'
#' @return The set of XGBoost relative importances between each transcription factor and the target gene.
#'
#' @seealso \code{\link{cor}}, \code{\link{XGBoostSolver}}
#'
#' @family solver methods
#'
#' @examples
#' # Load included Alzheimer's data, create a TReNA object with Bayes Spike as solver, and solve
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' XGBoost.solver <- XGBoostSolver(mtx.sub, target.gene, tfs)
#' tbl <- run(XGBoost.solver)

setMethod("run", "XGBoostSolver",

          function (obj){

              mtx <- getAssayData(obj)
              target.gene <- getTarget(obj)
              tfs <- getRegulators(obj)

              # Check that target gene and tfs are all part of the matrix
              stopifnot(target.gene %in% rownames(mtx))
              stopifnot(all(tfs %in% rownames(mtx)))

              # If given no tfs, return nothing
              if (length(tfs)==0) return(data.frame())

              # Don't handle tf self-regulation, so take target gene out of tfs
              deleters <- grep(target.gene, tfs)
              if(length(deleters) > 0){
                  tfs <- tfs[-deleters]
              }
              # If target gene was the only tf, then return nothing
              if(length(tfs)==0) return(NULL)

              x = t(mtx[tfs,,drop=FALSE])
              y = as.vector(t(mtx[target.gene,])) # Make target gene levels into a vector

              # call XGBoost here
              bst <- xgboost(data = x, label = y, eta = 0.1, nrounds = 100, verbose=0)
              tbl.importance <- xgb.importance(feature_names = colnames(x), model = bst)
              plot.importance <- xgb.plot.importance(importance_matrix = tbl.importance, plot=FALSE)

              # Fashion a data.frame with Gain, Cover, and Frequency for each transcription factor
              tbl <- data.frame(row.names = tbl.importance$Feature, tbl.importance[,-1])

              return(tbl)
          })
#----------------------------------------------------------------------------------------------------
