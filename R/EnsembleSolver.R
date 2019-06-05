#' Class EnsembleSolver
#'
#' @import methods
#' @import utils
#' @include Solver.R
#'
#' @name EnsembleSolver-class
#'
.EnsembleSolver <- setClass("EnsembleSolver",
                            contains="Solver",
                            slots = c(solverNames = "character",
                                      geneCutoff = "numeric",
                                      alpha.lasso = "numeric",
                                      alpha.ridge = "numeric",
                                      lambda.lasso = "numeric",
                                      lambda.ridge = "numeric",
                                      lambda.sqrt = "numeric",
                                      nCores.sqrt = "numeric",
                                      nOrderings.bayes = "numeric"
                                      )
                            )
#----------------------------------------------------------------------------------------------------
setGeneric("getSolverNames", signature = "obj", function(obj) standardGeneric("getSolverNames"))
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using an ensemble of solvers
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param targetGene A designated target gene that should be part of the mtx.assay data
#' @param candidateRegulators The designated set of transcription factors that could be associated
#' with the target gene
#' @param solverNames A character vector of strings denoting
#' @param geneCutoff A fraction (0-1) of the supplied candidate regulators to be included in the
#' fetaures output by the solver (default = 0.1)
#' @param alpha.lasso A fraction (0-1) denoting the LASSO-Ridge balance of the `glmnet` solver used
#' by the LASSO method (default = 0.9)
#' @param alpha.ridge A fraction (0-1) denoting the LASSO-Ridge balance of the `glmnet` solver used
#' by the Ridge method (default = 0)
#' @param lambda.lasso The penalty parameter for LASSO, used to determine how strictly to penalize
#' the regression coefficients. If none is supplied, this will be determined via permutation
#' testing (default = NULL).
#' @param lambda.ridge The penalty parameter for Ridge, used to determine how strictly to penalize
#' the regression coefficients. If none is supplied, this will be determined via permutation
#' testing (default = NULL).
#' @param lambda.sqrt The penalty parameter for square root LASSO, used to determine how strictly
#' to penalize the regression coefficients. If none is supplied, this will be determined via
#' permutation testing (default = NULL).
#' @param nCores.sqrt An integer denoting the number of computational cores to devote to the
#' square root LASSO solver, which is the slowest of the solvers (default = 4)
#' @param nOrderings.bayes An integer denoting the number of random starts to use for the Bayes
#' Spike method (default = 10)
#' @param quiet A logical denoting whether or not the solver should print output
#'
#' @return A Solver class object with Ensemble as the solver
#'
#' @seealso  \code{\link{solve.Ensemble}}, \code{\link{getAssayData}}
#'
#' @family Solver class objects
#'
#' @export
#'
#' @examples
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' ensemble.solver <- EnsembleSolver(mtx.sub, target.gene, tfs)


EnsembleSolver <- function(mtx.assay=matrix(), targetGene, candidateRegulators,
                           solverNames = c("lasso",
                                           "lassopv",
                                           "pearson",
                                           "randomForest",
                                           "ridge",
                                           "spearman"),
                           geneCutoff = 0.1,
                           alpha.lasso = 0.9,
                           alpha.ridge = 0.0,
                           lambda.lasso = numeric(0),
                           lambda.ridge = numeric(0),
                           lambda.sqrt = numeric(0),
                           nCores.sqrt = 4,
                           nOrderings.bayes = 10,
                           quiet=TRUE)
{
    if(any(grepl(targetGene, candidateRegulators)))
        candidateRegulators <- candidateRegulators[-grep(targetGene, candidateRegulators)]

    candidateRegulators <- intersect(candidateRegulators, rownames(mtx.assay))
    stopifnot(length(candidateRegulators) > 0)

    # Send a warning if there's a row of zeros
    if(!is.na(max(mtx.assay)) & any(rowSums(mtx.assay) == 0))
        warning("One or more gene has zero expression; this may cause difficulty when using Bayes Spike. You may want to try 'lasso' or 'ridge' instead.")

    obj <- .EnsembleSolver(Solver(mtx.assay = mtx.assay,
                                  targetGene = targetGene,
                                  candidateRegulators = candidateRegulators,
                                  quiet = quiet),
                           solverNames = solverNames,
                           geneCutoff = geneCutoff,
                           alpha.lasso = alpha.lasso,
                           alpha.ridge = alpha.ridge,
                           lambda.lasso = lambda.lasso,
                           lambda.ridge = lambda.ridge,
                           lambda.sqrt = lambda.sqrt,
                           nCores.sqrt = nCores.sqrt,
                           nOrderings.bayes = nOrderings.bayes)
    obj

} # EnsembleSolver, the constructor
#----------------------------------------------------------------------------------------------------
#' Show the Ensemble Solver
#'
#' @rdname show.EnsembleSolver
#' @aliases show.EnsembleSolver
#'
#' @param object An object of the class EnsembleSolver
#'
#' @return A truncated view of the supplied object
#'
#' @examples
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' ensemble.solver <- EnsembleSolver(mtx.sub, target.gene, tfs)
#' show(ensemble.solver)

setMethod('show', 'EnsembleSolver',

          function(object) {
              regulator.count <- length(getRegulators(object))
              if(regulator.count > 10){
                  regulatorString <- paste(getRegulators(object)[1:10], collapse=",")
                  regulatorString <- sprintf("%s...", regulatorString);
              }
              else
                  regulatorString <- paste(getRegulators(object), collapse=",")

              msg = sprintf("EnsembleSolver with mtx.assay (%d, %d), targetGene %s, %d candidate regulators %s,  solvers: %s",
                            nrow(getAssayData(object)), ncol(getAssayData(object)),
                            getTarget(object), regulator.count, regulatorString,
                            paste(object@solverNames,collapse = ", "))
              cat (msg, '\n', sep='')
          })
#----------------------------------------------------------------------------------------------------
#' Retrieve the solver names from an EnsembleSolver object
#'
#' @rdname getSolverNames
#' @aliases getSolverNames
#'
#' @param obj An object of class Solver
#'
#' @return The vector of solver names associated with an EnsembleSolver object
#'
#' @examples
#' # Create a Solver object using the included Alzheimer's data and retrieve the regulators
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' targetGene <- "MEF2C"
#' candidateRegulators <- setdiff(rownames(mtx.sub), targetGene)
#' solver <- EnsembleSolver(mtx.sub, targetGene, candidateRegulators,
#' solverNames = c("lasso","randomForest"))
#' solver.names <- getSolverNames(solver)

#' @export

setMethod("getSolverNames", "EnsembleSolver",

          function (obj){
              obj@solverNames
          })
#----------------------------------------------------------------------------------------------------
#' Run the Ensemble Solver
#'
#' @name run,EnsembleSolver-method
#' @rdname solve.Ensemble
#' @aliases run.EnsembleSolver solve.Ensemble
#'
#' @description Given a TReNA object with Ensemble as the solver and a list of solvers
#' (default = "default.solvers"), estimate coefficients for each transcription factor
#' as a predictor of the target gene's expression level. The final scores for the ensemble
#' method combine all specified solvers to create a composite score for each transcription factor.
#' This method should be called using the \code{\link{solve}} method on an appropriate TReNA object.
#'
#' @param obj An object of class Solver with "ensemble" as the solver string
#'
#' @return A data frame containing the scores for all solvers and two composite scores
#' relating the target gene to each transcription factor. The two new scores are:
#' @details
#' \itemize{
#'   \item{concordance}{a composite score}
#'   \item{pcaMax}{a composite of the principal components of the individual solver scores}
#' }
#'
#' @seealso \code{\link{EnsembleSolver}}
#'
#' @family solver methods
#'
#' @examples
#' \dontrun{
#' # Load included Alzheimer's data, create an Ensemble object with default solvers, and solve
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)[1:30]
#' ensemble.solver <- EnsembleSolver(mtx.sub, target.gene, tfs)
#' tbl <- run(ensemble.solver)
#'
#' # Solve the same problem, but supply extra arguments that change alpha for LASSO to 0.8 and also
#' # Change the gene cutoff from 10% to 20%
#' ensemble.solver <- EnsembleSolver(mtx.sub, target.gene, tfs, geneCutoff = 0.2, alpha.lasso = 0.8)
#' tbl <- run(ensemble.solver)
#'
#' # Solve the original problem with default cutoff and solver parameters, but use only 4 solvers
#' ensemble.solver <- EnsembleSolver(mtx.sub, target.gene, tfs,
#' solverNames = c("lasso", "pearson", "ridge"))
#' tbl <- run(ensemble.solver)
#' }
#'
#' @export

setMethod("run", "EnsembleSolver",

          function(obj){

              mtx <- getAssayData(obj)
              target.gene <- getTarget(obj)
              tfs <- getRegulators(obj)
              gene.cutoff <- obj@geneCutoff

              # Create a list of solvers and a list for solutions
              out.list <- list()
              solver.list <- tolower(getSolverNames(obj))

              # Intersect with the accepted solvers
              accepted.solvers <- c("bayesspike", "lassopv", "lasso", "pearson",
                                    "ridge", "randomforest", "spearman", "sqrtlasso")
              not.accepted <- setdiff(solver.list, accepted.solvers)
              solver.list <- intersect(solver.list, accepted.solvers)

              # If there's no valid solvers, exit with an error
              if(length(solver.list) == 0){
                  stop("No valid solvers supplied;
                       Run getAvailableSolvers() to see all available solvers")
              }

              # If there's any invalid solvers, throw a warning
              if(length(not.accepted) > 0){
                  not.accepted <- paste(not.accepted,collapse = ", ")
                  warn.msg <- sprintf("Invalid solvers (%s) supplied, running with remaining solver(s);
                                       Run getAvailableSolvers() to see all available solvers",
                                      not.accepted)
                  warning(warn.msg)
              }

              # If there's only 1 solver, catch it, send a warning, and run THAT solver
              if(length(solver.list) == 1){

                  # Capture the correct solver
                  solver <- switch(solver.list,
                                   "lasso" = LassoSolver(mtx, target.gene, tfs,
                                                         alpha = obj@alpha.lasso, lambda = obj@lambda.lasso),
                                   "randomforest" = RandomForestSolver(mtx, target.gene, tfs),
                                   "bayesspike" = BayesSpikeSolver(mtx, target.gene, tfs,
                                                                   nOrderings = obj@nOrderings.bayes),
                                   "pearson" = PearsonSolver(mtx, target.gene, tfs),
                                   "spearman" = SpearmanSolver(mtx, target.gene, tfs),
                                   "sqrtlasso" = SqrtLassoSolver(mtx, target.gene, tfs,
                                                                 lambda = obj@lambda.sqrt, nCores = obj@nCores.sqrt),
                                   "lassopv" = LassoPVSolver(mtx, target.gene, tfs),
                                   "ridge" = RidgeSolver(mtx, target.gene, tfs,
                                                         alpha = obj@alpha.ridge, lambda = obj@lambda.ridge))

                  # Send a specific warning message
                  warn.message <- sprintf("Only one solver(%s) was provided. Running %s instead",
                                          solver.list, class(solver)[1])
                  warning(warn.message)

                  # Run the solver
                  tbl.single <- run(solver)

                  # Return the output of the solver
                  return(tbl.single)
              }


              for(i in 1:length(solver.list)){
                  # Find the correct solver
                  solver <- switch(solver.list[i],
                                   "lasso" = LassoSolver(mtx, target.gene, tfs,
                                                         alpha = obj@alpha.lasso, lambda = obj@lambda.lasso),
                                   "randomforest" = RandomForestSolver(mtx, target.gene, tfs),
                                   "bayesspike" = BayesSpikeSolver(mtx, target.gene, tfs,
                                                                   nOrderings = obj@nOrderings.bayes),
                                   "pearson" = PearsonSolver(mtx, target.gene, tfs),
                                   "spearman" = SpearmanSolver(mtx, target.gene, tfs),
                                   "sqrtlasso" = SqrtLassoSolver(mtx, target.gene, tfs,
                                                                 lambda = obj@lambda.sqrt, nCores = obj@nCores.sqrt),
                                   "lassopv" = LassoPVSolver(mtx, target.gene, tfs),
                                   "ridge" = RidgeSolver(mtx, target.gene, tfs,
                                                         alpha = obj@alpha.ridge, lambda = obj@lambda.ridge))

                  # Solve each Solver object and save it to the output list
                  out.list[[i]] <- run(solver)
                  names(out.list)[i] <- paste("out",tolower(solver.list[[i]]),sep=".")
              }

              # Output lasso with beta
              if("lasso" %in% tolower(solver.list)){
                  out.list$out.lasso$gene <- rownames(out.list$out.lasso)
                  out.list$out.lasso <- out.list$out.lasso[, c("beta","gene")]
                  rownames(out.list$out.lasso) <- NULL
                  names(out.list$out.lasso) <- c("betaLasso", "gene")
                  lasso.med <- stats::median(out.list$out.lasso$betaLasso)
                  lasso.scale <- stats::mad(out.list$out.lasso$betaLasso)
              }

              # Output randomforest IncNodePurity
              if("randomforest" %in% tolower(solver.list)){
                  out.list$out.randomforest$gene <- rownames(out.list$out.randomforest)
                  out.list$out.randomforest <- out.list$out.randomforest[, c("IncNodePurity","gene")]
                  rownames(out.list$out.randomforest) <- NULL
                  names(out.list$out.randomforest) <- c("rfScore", "gene")
                  randomforest.med <- stats::median(out.list$out.randomforest$rfScore)
                  randomforest.scale <- sqrt(mean(
                      out.list$out.randomforest$rfScore*out.list$out.randomforest$rfScore))
              }

              # Output the z-score from bayesspike
              if("bayesspike" %in% tolower(solver.list)){
                  out.list$out.bayesspike$gene <- rownames(out.list$out.bayesspike)
                  rownames(out.list$out.bayesspike) <- NULL
                  out.list$out.bayesspike <- out.list$out.bayesspike[, c("z", "gene")]
                  names(out.list$out.bayesspike) <- c("bayesScore", "gene")
                  bayesspike.med <- stats::median(out.list$out.bayesspike$bayesScore)
                  bayesspike.scale <- stats::mad(out.list$out.bayesspike$bayesScore)
              }

              # Pearson
              if("pearson" %in% tolower(solver.list)){
                  out.list$out.pearson$gene <- rownames(out.list$out.pearson)
                  rownames(out.list$out.pearson) <- NULL
                  names(out.list$out.pearson) <- c("pearsonCoeff","gene")
                  pearson.med <- stats::median(out.list$out.pearson$pearsonCoeff)
                  pearson.scale <- stats::mad(out.list$out.pearson$pearsonCoeff)
              }

              #Spearman
              if("spearman" %in% tolower(solver.list)){
                  out.list$out.spearman$gene <- rownames(out.list$out.spearman)
                  rownames(out.list$out.spearman) <- NULL
                  names(out.list$out.spearman) <- c("spearmanCoeff", "gene")
                  spearman.med <- stats::median(out.list$out.spearman$spearmanCoeff)
                  spearman.scale <- stats::mad(out.list$out.spearman$spearmanCoeff)
              }

              #LassoPV
              if("lassopv" %in% tolower(solver.list)){
                  out.list$out.lassopv$gene <- rownames(out.list$out.lassopv)
                  rownames(out.list$out.lassopv) <- NULL
                  out.list$out.lassopv <- out.list$out.lassopv[, c("pValues","gene")]
                  names(out.list$out.lassopv) <- c("lassoPValue", "gene")
                  p.log10 <- -log10(out.list$out.lassopv$lassoPValue)
                  lassopv.med <- stats::median(p.log10)
                  lassopv.scale <- sqrt(mean(p.log10*p.log10))
              }

              #SqrtLasso
              if("sqrtlasso" %in% tolower(solver.list)){
                  out.list$out.sqrtlasso$gene <- rownames(out.list$out.sqrtlasso)
                  rownames(out.list$out.sqrtlasso) <- NULL
                  out.list$out.sqrtlasso <- out.list$out.sqrtlasso[, c("beta", "gene")]
                  names(out.list$out.sqrtlasso) <- c("betaSqrtLasso", "gene")
                  sqrtlasso.med <- stats::median(out.list$out.sqrtlasso$betaSqrtLasso)
                  sqrtlasso.scale <- stats::mad(out.list$out.sqrtlasso$betaSqrtLasso)
              }

              #Ridge
              if("ridge" %in% tolower(solver.list)){
                  out.list$out.ridge$gene <- rownames(out.list$out.ridge)
                  out.list$out.ridge <- out.list$out.ridge[, c("beta","gene")]
                  rownames(out.list$out.ridge) <- NULL
                  names(out.list$out.ridge) <- c("betaRidge", "gene")
                  ridge.med <- stats::median(out.list$out.ridge$betaRidge)
                  ridge.scale <- stats::mad(out.list$out.ridge$betaRidge)
              }

              # Grab all genes for each solver to start with
              how.many <- length(tfs)
              top.list <- lapply(out.list, function(x) head(x$gene, how.many))
              all.genes <- unique(as.character(unlist(top.list)))

              # Run same thing in a while loop until the cutoff or 10 is reached
              while(length(all.genes) > gene.cutoff * length(tfs) & length(all.genes) > 10){

                  how.many <- round(0.75*how.many)
                  top.list <- lapply(out.list, function(x) utils::head(x$gene, how.many))
                  all.genes <- unique(as.character(unlist(top.list)))
              }

              # Pull out the specified genes
              sub.list <- list()
              for(i in 1:length(out.list)){

                  sub.list[[i]] <- subset(out.list[[i]], out.list[[i]]$gene %in% all.genes)
              }

              # Merge the tables
              tbl.all <- merge(sub.list[[1]], sub.list[[2]], by = "gene", all = TRUE)
              if(length(sub.list) > 2){
                  for(i in 3:(length(sub.list))){
                      tbl.all <- merge(tbl.all, sub.list[[i]], by = "gene", all = TRUE)
                      }
                  }

              tbl.all[is.na(tbl.all)] <- 0

#              # Replace missing values and scale the data
#              # Use the *.med and *.scale values to center/scale everything
#              tbl.scale <- tbl.all[,-1]
#
#              if("lassoPValue" %in% names(tbl.scale)){
#                  tbl.scale$lassoPValue <- -log10(tbl.scale$lassoPValue)
#                  tbl.scale$lassoPValue <- scale(tbl.scale$lassoPValue,
#                                                   center = lassopv.med,
#                                                   scale = lassopv.scale)
#              }
#
#              if("betaLasso" %in% names(tbl.scale)){
#                  tbl.scale$betaLasso <- scale(tbl.scale$betaLasso,
#                                                center = lasso.med,
#                                                scale = lasso.scale)
#              }
#
#              if("betaRidge" %in% names(tbl.scale)){
#                  tbl.scale$betaRidge <- scale(tbl.scale$betaRidge,
#                                                center = ridge.med,
#                                                scale = ridge.scale)
#              }
#
#              if("pearsonCoeff" %in% names(tbl.scale)){
#                  tbl.scale$pearsonCoeff <- scale(tbl.scale$pearsonCoeff,
#                                                   center = pearson.med,
#                                                   scale = pearson.scale)
#              }
#
#              if("spearmanCoeff" %in% names(tbl.scale)){
#                  tbl.scale$spearmanCoeff <- scale(tbl.scale$spearmanCoeff,
#                                                    center = spearman.med,
#                                                    scale = spearman.scale)
#              }
#
#              if("betaSqrtLasso" %in% names(tbl.scale)){
#                  tbl.scale$betaSqrtLasso <- scale(tbl.scale$betaSqrtLasso,
#                                                    center = sqrtlasso.med,
#                                                    scale = sqrtlasso.scale)
#              }
#
#              if("bayesScore" %in% names(tbl.scale)){
#                  tbl.scale$bayesScore <- scale(tbl.scale$bayesScore,
#                                             center = bayesspike.med,
#                                             scale = bayesspike.scale)
#              }
#
#              if("rfScore" %in% names(tbl.scale)){
#                  tbl.scale$rfScore <- scale(tbl.scale$rfScore,
#                                              center = randomforest.med,
#                                              scale = randomforest.scale)
#              }
#
#              rownames(tbl.scale) <- tbl.all$gene
#
#              tbl.augmented <- try(.addEnsembleScores(tbl.scale, tbl.all), silent = TRUE)
#
#              # If you get new scores, add them;
#              # Else, just keep the old table and throw a warning
#
#              if(class(tbl.augmented) == "try-error"){
#                  warning("The signal strength of ensemble of solvers is too weak to support
# composite scores ('pcaMax' and 'concordance' in the model output table. This is a classic
# "large n, small m" problem that could be rectified by providing more samples")
#                  tbl.all$pcaMax <- NA
#                  tbl.all$concordance <- NA
#              } else {
#                  tbl.all <- tbl.augmented
#              }
#
#              # Regardless of output, return the table of scores

          return(tbl.all)
          })
#----------------------------------------------------------------------------------------------------
#.addEnsembleScores <- function(tbl.scale, tbl.all) {
#
#    # Compute the scaled "concordance score"
#    pca <- stats::prcomp(tbl.scale, center=FALSE, scale.=FALSE)
#
#    pca$x <- pca$x / sqrt(length(which(pca$sdev > 0.1)))
#    concordance <- apply(pca$x[, pca$sdev > 0.1, drop=FALSE], 1,
#                         function(x) {sqrt(mean((2*atan(x)/pi)^2))})
#    concordance <- as.data.frame(concordance)
#    concordance$gene <- rownames(concordance)
#    rownames(concordance) <- NULL
#    tbl.all <- merge(tbl.all, concordance, by = "gene", all = TRUE)
#
#    # Transform via PCA and compute the pcaMax score
#    pcaMax <- apply(pca$x[, pca$sdev > 0.1, drop=FALSE],1, function(x) {sqrt(mean(x*x))})
#    pcaMax <- as.data.frame(pcaMax)
#    pcaMax$gene <- rownames(pcaMax)
#    rownames(pcaMax) <- NULL
#    tbl.all <- merge(tbl.all, pcaMax, by = "gene", all = TRUE)
#
#    # Sort by pcaMax
#    tbl.all <- tbl.all[order(tbl.all$pcaMax, decreasing = TRUE),]
#
#    return(tbl.all)
#}
#----------------------------------------------------------------------------------------------------
