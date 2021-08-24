#----------------------------------------------------------------------------------------------------
#' An S4 class to represent a Bicor solver
#'
#' @include Solver.R
#' @import methods
#' @importFrom WGCNA bicor
#'
#' @name BicorSolver-class
#'

.BicorSolver <- setClass ("BicorSolver",contains = "Solver")
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using  Bicor correlation coefficients as the solver
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param targetGene A designated target gene that should be part of the mtx.assay data
#' @param candidateRegulators The designated set of transcription factors that could be associated
#' with the target gene
#' @param quiet A logical denoting whether or not the solver should print output
#'
#' @return A Solver class object with Bicor correlation coefficients as the solver
#'
#' @seealso  \code{\link{solve.Bicor}}, \code{\link{getAssayData}}
#'
#' @family Solver class objects
#'
#' @export
#'
#' @examples
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' bicor.solver <- BicorSolver(mtx.sub, target.gene, tfs)

BicorSolver <- function(mtx.assay = matrix(), targetGene, candidateRegulators, quiet=TRUE)
{
    # Remove the targetGene from candidateRegulators
    if(any(grepl(targetGene, candidateRegulators)))
        candidateRegulators <- candidateRegulators[-grep(targetGene, candidateRegulators)]

    # Check to make sure the matrix contains some of the candidates
    candidateRegulators <- intersect(candidateRegulators, rownames(mtx.assay))
    stopifnot(length(candidateRegulators) > 0)

    obj <- .BicorSolver(Solver(mtx.assay=mtx.assay,
                                 quiet=quiet,
                                 targetGene=targetGene,
                                 candidateRegulators=candidateRegulators))

    # Send a warning if there's a row of zeros
    if(!is.na(max(mtx.assay)) & any(rowSums(mtx.assay) == 0))
        warning("One or more gene has zero expression; this may yield 'NA' results and warnings when using Bicor correlations")

    obj

} #BicorSolver, the constructor
#----------------------------------------------------------------------------------------------------
#' Show the Bicor Solver
#'
#' @rdname show.BicorSolver
#' @aliases show.BicorSolver
#'
#' @param object An object of the class BicorSolver
#'
#' @return A truncated view of the supplied object
#'
#' @examples
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' bicor.solver <- BicorSolver(mtx.sub, target.gene, tfs)
#' show(bicor.solver)

setMethod('show', 'BicorSolver',

          function(object) {
              regulator.count <- length(getRegulators(object))
              if(regulator.count > 10){
                  regulatorString <- paste(getRegulators(object)[1:10], collapse=",")
                  regulatorString <- sprintf("%s...", regulatorString);
              }
              else
                  regulatorString <- paste(getRegulators(object), collapse=",")

              msg = sprintf("BicorSolver with mtx.assay (%d, %d), targetGene %s, %d candidate regulators %s",
                            nrow(getAssayData(object)), ncol(getAssayData(object)),
                            getTarget(object), regulator.count, regulatorString)
              cat (msg, '\n', sep='')
          })
#----------------------------------------------------------------------------------------------------
#' Run the Bicor Solver
#'
#' @rdname solve.Bicor
#' @aliases run.BicorSolver solve.Bicor
#'
#' @description Given a BicorSolver object, use the \code{\link{cor}} function
#' to estimate coefficients for each transcription factor as a predictor of the target gene's
#' expression level.
#'
#' @param obj An object of class BicorSolver
#'
#' @return The set of Bicor Correlation Coefficients between each transcription factor and the target gene.
#'
#' @seealso \code{\link{cor}}, \code{\link{BicorSolver}}
#'
#' @family solver methods
#'
#' @examples
#' # Load included Alzheimer's data, create a TReNA object with Bayes Spike as solver, and solve
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' bicor.solver <- BicorSolver(mtx.sub, target.gene, tfs)
#' tbl <- run(bicor.solver)

setMethod("run", "BicorSolver",

          function (obj){

              mtx <- getAssayData(obj)
              target.gene <- getTarget(obj)
              tfs <- getRegulators(obj)

              # Check that target gene and tfs are all part of the matrix
              stopifnot(target.gene %in% rownames(mtx))
              stopifnot(all(tfs %in% rownames(mtx)))
              # If given no tfs, return nothing
              if (length(tfs)==0) return(NULL)

              # Don't handle tf self-regulation, so take target gene out of tfs
              deleters <- grep(target.gene, tfs)
              if(length(deleters) > 0){
                  tfs <- tfs[-deleters]
              }
              # If target gene was the only tf, then return nothing
              if(length(tfs)==0) return(NULL)

              x = t(mtx[tfs,,drop=FALSE])
              y = as.vector(t(mtx[target.gene,])) # Make target gene levels into a vector

              # Calculate Bicor correlation coefficients
              fit <- WGCNA::bicor(x=x, y=y, use="pairwise.complete.obs")

              # Return the coefficients as a data frame
              tbl <- data.frame(row.names = rownames(fit)[order(abs(fit), decreasing = TRUE)],
                                coefficient = fit[order(abs(fit), decreasing = TRUE)])

              return(tbl)
          })
#----------------------------------------------------------------------------------------------------
