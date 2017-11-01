#----------------------------------------------------------------------------------------------------
#' An S4 class to represent a Bayes Spike solver
#'
#' @import vbsr
#' @import methods
#' 
#' @include Solver.R
#' @name BayesSpikeSolver-class

.BayesSpikeSolver <- setClass("BayesSpikeSolver",
                              contains="Solver",
                              slots = c(nOrderings = "numeric")
                              )
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using the Bayes Spike Solver
#' 
#' @param mtx.assay An assay matrix of gene expression data
#' @param targetGene A designated target gene that should be part of the mtx.assay data
#' @param candidateRegulators The designated set of transcription factors that could be associated
#' with the target gene
#' @param nOrderings An integer denoting the number of random starts to use for the Bayes Spike
#' method (default = 10)
#' @param quiet A logical denoting whether or not the solver should print output
#'
#' @export
#' 
#' @return A Solver class object with Bayes Spike as the solver
#'
#' @family Solver class objects
#' 
#' @seealso  \code{\link{solve.BayesSpike}}, \code{\link{getAssayData}}
#'
#' @examples
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' bayes.solver <- BayesSpikeSolver(mtx.sub, target.gene, tfs)

BayesSpikeSolver <- function(mtx.assay=matrix(), targetGene, candidateRegulators,
                             nOrderings = 10, quiet=TRUE)
{
    if(any(grepl(targetGene, candidateRegulators)))        
        candidateRegulators <- candidateRegulators[-grep(targetGene, candidateRegulators)]    
    
    candidateRegulators <- intersect(candidateRegulators, rownames(mtx.assay))    
    stopifnot(length(candidateRegulators) > 0)    
    
    obj <- .BayesSpikeSolver(Solver(mtx.assay=mtx.assay,                             
                                    targetGene = targetGene,
                                    candidateRegulators = candidateRegulators,
                                    quiet = quiet),
                             nOrderings = nOrderings)
    
    # Send a warning if there's a row of zeros
    if(!is.na(max(mtx.assay)) & any(rowSums(mtx.assay) == 0))
        warning("One or more gene has zero expression; this may cause difficulty when using Bayes Spike. You may want to try 'lasso' or 'ridge' instead.")
    
    obj
    
} # BayesSpikeSolver, the constructor
#----------------------------------------------------------------------------------------------------
#' Show the Bayes Spike Solver
#' 
#' @rdname show.BayesSpikeSolver
#' @aliases show.BayesSpikeSolver
#'
#' @param object An object of the class BayesSpikeSolver
#'
#' @return A truncated view of the supplied object
#'
#' @examples
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' bayes.solver <- BayesSpikeSolver(mtx.sub, target.gene, tfs)
#' show(bayes.solver)

setMethod('show', 'BayesSpikeSolver',
          
          function(object) {
              regulator.count <- length(getRegulators(object))
              if(regulator.count > 10){
                  regulatorString <- paste(getRegulators(object)[1:10], collapse=",")
                  regulatorString <- sprintf("%s...", regulatorString);
              }
              else
                  regulatorString <- paste(getRegulators(object), collapse=",")
              
              msg = sprintf("BayesSpikeSolver with mtx.assay (%d, %d), targetGene %s, %d candidate regulators %s, %d orderings",
                            nrow(getAssayData(object)), ncol(getAssayData(object)),
                            getTarget(object), regulator.count, regulatorString, object@nOrderings)
              cat (msg, '\n', sep='')
          })
#----------------------------------------------------------------------------------------------------
#' Run the Bayes Spike Solver
#'
#' @rdname solve.BayesSpike
#' @aliases run.BayesSpikeSolver solve.BayesSpike
#' 
#' @description Given a TReNA object with Bayes Spike as the solver, use the \code{\link{vbsr}}
#' function to estimate coefficients for each transcription factor as a predictor of the target
#' gene's expression level.
#'
#' @param obj An object of the class BayesSpikeSolver
#' 
#' @return A data frame containing the coefficients relating the target gene to each transcription factor,
#' plus other fit parameters
#'
#' @seealso \code{\link{vbsr}}, \code{\link{BayesSpikeSolver}}
#'
#' @family solver methods
#' 
#' @examples
#' \dontrun{
#' # Load included Alzheimer's data, create a TReNA object with Bayes Spike as solver, and solve
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' bayes.solver <- BayesSpikeSolver(mtx.sub, target.gene, tfs)
#' tbl <- run(bayes.solver)
#'
#' # Solve the same Alzheimer's problem, but this time set the number of random starts to 100
#' bayes.solver <- BayesSpikeSolver(mtx.sub, target.gene, tfs, nOrderings = 100)
#' tbl <- run(bayes.solver)
#' }

setMethod("run", "BayesSpikeSolver",

          function (obj){
              
              mtx <- getAssayData(obj)    
              target.gene <- getTarget(obj)
              tfs <- getRegulators(obj)
              stopifnot(target.gene %in% rownames(mtx))      
              stopifnot(all(tfs %in% rownames(mtx)))
              
              # we don't try to handle tf self-regulation      
              deleters <- grep(target.gene, tfs)      
              if(length(deleters) > 0){          
                  tfs <- tfs[-deleters]          
                  message(sprintf("BayesSpikeSolver removing target.gene from candidate regulators: %s", target.gene))          
              }      
              
              features <- t(mtx[tfs, ])      
              target <- as.numeric(mtx[target.gene,])      
              result <- vbsr(target, features, family='normal', n_orderings = obj@nOrderings)
              
              # Add Pearson coefficient and add a "score"
              tbl.out <- data.frame(beta=result$beta, pval=result$pval, z=result$z, post=result$post)      
              rownames(tbl.out) <- tfs      
              tbl.out$score <- -log10(tbl.out$pval)      
              tbl.out <- tbl.out[order(tbl.out$score, decreasing=TRUE),]
              tbl.out
              
          })
#----------------------------------------------------------------------------------------------------
