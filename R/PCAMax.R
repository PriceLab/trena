#' Class PCAMax
#'
#' @import methods
#' @import utils
#'
#' @name PCAMax
#'
.pcaMax <- setClass("PCAMax",
                    slots=c(tbl="data.frame",
                            solverNames="character",
                            state="environment")
                    )

#------------------------------------------------------------------------------------------------------------------------
setGeneric("normalizeModel",   signature="obj", function(obj, normalizing.max=10) standardGeneric ("normalizeModel"))
setGeneric("addStats",         signature="obj", function(obj) standardGeneric ("addStats"))
setGeneric("getCoverage",      signature="obj", function(obj) standardGeneric ("getCoverage"))
#------------------------------------------------------------------------------------------------------------------------
#' Create a PCAMax object from a data.frame produced by the EnsembleSolver
#'
#' @param tbl a data.frame
#'
#' @return a PCAMax object
#'
#' @seealso  \code{\link{EnsembleSolver}}
#'
#' @export
#'
PCAMax <- function(tbl, tfIdentifierColumnName="tf.hgnc")
{
    stopifnot(tfIdentifierColumnName %in% colnames(tbl))
    standardSolverNames <- c("betaLasso", "lassoPValue", "pearsonCoeff", "rfScore",
                             "betaRidge", "spearmanCoeff", "betaSqrtLasso")

    coi <- intersect(standardSolverNames, colnames(tbl))
    if(length(coi) < 3)
      stop(sprintf("insufficient dimensions,  %d in EnsembleSolver data.frame to calculate pcaMax", length(coi)))

    tbl.trimmed <- tbl[, coi]
    rownames(tbl.trimmed) <- tbl[, tfIdentifierColumnName]
    state <- new.env(parent=emptyenv())
    state$coverage <- 0
    .pcaMax(tbl=tbl.trimmed, solverNames=standardSolverNames, state=state)

} # PCAMax ctor
#------------------------------------------------------------------------------------------------------------------------
#' transform a specific column to fit normal distribution
#'
#' @rdname normalizeModel
#' @aliases normalizeModel
#'
#' @param obj An object of the class PCAMax
#' @param max numeric, a maximum value for the normalized distrubtions
#'
#' @return a normalized matrix, each column treated separately
#'
#' @export
#'
setMethod('normalizeModel', 'PCAMax',

     function(obj, normalizing.max=10){

        mtx <- as.matrix(obj@tbl)
        for(colname in colnames(mtx)){
           vec <- as.numeric(mtx[, colname])
           if(colname %in% c("betaLasso", "betaRidge", "betaSqrtLasso"))
              vec.norm <- .normalize_betaValues(vec, normalizing.max)
           if(colname == "lassoPValue")
              vec.norm <- .normalize_pval(vec, normalizing.max)
           if(colname == "rfScore")
              vec.norm <- .normalize_randomForest(vec, normalizing.max)
           if(colname %in% c("spearmanCoeff", "pearsonCoeff"))
              vec.norm <- .normalize_correlationValues(vec, normalizing.max)
           mtx[, colname] <- vec.norm
           } # for colname
        obj@state$normalizedMatrix <- mtx
        mtx
        })

#------------------------------------------------------------------------------------------------------------------------
#' add PCA-based summary stats on all TFs in the model
#'
#' @rdname addStats
#' @aliases addStats
#'
#' @param obj An object of the class PCAMax
#'
#' @return  the original model with extra columns: pcaMax, cov, PC1, PC2,
#'
#' @export
#'
setMethod('addStats', 'PCAMax',

     function(obj){
        mtx <- obj@state$normalizedMatrix
        pca <- prcomp(mtx)
        mtx.pca <- pca$x
         # coverage: the variance covered by the first two components
        obj@state$coverage <- sum(pca$sdev[1:2]) / sum(pca$sdev)
        mean <- apply(mtx, 1, function(row) mean(abs(row)))
        mtx.2 <- cbind(mtx, mean)
        rms <- apply(mtx, 1, function(row) sqrt(mean(row * row)))
        pcaMax <- as.numeric(apply(mtx.pca, 1, function(row) sqrt((row[1] * row[1]) + (row[2] * row[2]))))
        mtx.3 <- cbind(mtx.2, rms)
        mtx.4 <- cbind(mtx.3, mtx.pca[, 1:2])
        mtx.5 <- cbind(mtx.4, pcaMax)
          #  use coefficient of variation to desribe the dispersion (agreement) of the data
        covar <- apply(mtx, 1, function(row) sd(abs(row))/mean(abs(row)))
        mtx.6 <- cbind(mtx.5, covar)
        mtx.7 <- mtx.6[order(mtx.6[, "pcaMax"], decreasing=TRUE),]
        obj@state$normalizedMatrixWithStats <- mtx.7
        mtx.7
        })

#------------------------------------------------------------------------------------------------------------------------
#' what percentage of the variance is captured in the first two principal components?
#'
#' @rdname getCoverage
#' @aliases getCoverage
#'
#' @param obj An object of the class PCAMax
#'
#' @return a number between zero and one
#'
#' @export
#'
setMethod('getCoverage', 'PCAMax',

     function(obj){
        obj@state$coverage
     })

#------------------------------------------------------------------------------------------------------------------------
.normalize_randomForest <- function(vector, normalizing.max=10)
{
   third.quartile.normalizing.max <- 0.75 * normalizing.max
   third.quartile <- fivenum(vector)[4]
   actual.max <- third.quartile
   scale <- third.quartile.normalizing.max/actual.max
   vector * scale

} # .normalize_randomForest
#------------------------------------------------------------------------------------------------------------------------
.normalize_pval <- function(vector, normalizing.max=10)
{
   stopifnot(all(vector >= 0))
   x <- -log10(vector +  .Machine$double.xmin)  #   ensure no zero values
   actual.max <- max(x)
   scale <- normalizing.max/actual.max
   x * scale

} # .normalize_pval
#------------------------------------------------------------------------------------------------------------------------
.normalize_betaValues <- function(vector, normalizing.max=10)
{
   actual.max <- max(abs(vector))
   scale <- normalizing.max/actual.max
   vec.new <- vector * scale
   vec.new

} # .normalize_betaValues
#------------------------------------------------------------------------------------------------------------------------
.normalize_correlationValues <- function(vector, normalizing.max=10)
{
   actual.max <- 1
   scale <- normalizing.max/actual.max
   vector * scale

} # .normalize_correlationValues
#------------------------------------------------------------------------------------------------------------------------
#              # Replace missing values and scale the data
#              # Use the *.med and *.scale values to center/scale everything
#              tbl.all[is.na(tbl.all)] <- 0
#              tbl.scale <- tbl.all[-1]
#
#              if(
# "lassoPValue" %in% names(tbl.scale)){
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
##------------------------------------------------------------------------------------------------------------------------
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
##------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------

