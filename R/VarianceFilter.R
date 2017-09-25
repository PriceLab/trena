#' @title Create a VarianceFilter object
#'
#' @description
#' A VarianceFilter object allows for filtering based on the variance of a target gene in relation to
#' other genes in the assay matrix. Using its associated \code{getCandidates} method, a VarianceFilter
#' object can be used to filter a list of possible transcription factors to those within a given
#' range of the variance of a supplied target gene.
#'
#' @include CandidateFilter.R
#' @import methods
#'
#' @name VarianceFilter-class
#' @rdname VarianceFilter-class
#' @aliases VarianceFilter

#----------------------------------------------------------------------------------------------------
.VarianceFilter <- setClass("VarianceFilter",
                            contains = "CandidateFilter",
                            slots = c(mtx.assay = "matrix",
                                      targetGene = "character",
                                      varSize = "numeric")
                            )

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
#' @rdname VarianceFilter-class
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param varSize A user-specified fraction (0-1) of the target gene variance to use as a filter
#' @param targetGene A designated target gene that must be part of the mtx.assay data
#' @param quiet A logical denoting whether or not the solver should print output
#'
#' @return A CandidateFilter class object with variance as the filtering method
#'
#' @export
#' 
#' @seealso \code{\link{getCandidates-VarianceFilter}}
#'
#' @return An object of the VarianceFilter class
#' 
#' @family Filtering Objects
#'
#' @examples
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' variance.filter <- VarianceFilter(mtx.assay = mtx.sub, targetGene = "MEF2C")

VarianceFilter <- function(mtx.assay=matrix(), targetGene, varSize = 0.5, quiet=TRUE)
{
    .VarianceFilter(mtx.assay = mtx.assay,
                    targetGene = targetGene,
                    varSize = varSize,
                    quiet = quiet)

} # VarianceFilter, the constructor
#----------------------------------------------------------------------------------------------------
#' Get candidate genes using the variance filter
#'
#' @aliases getCandidates-VarianceFilter
#'
#' @param obj An object of class VarianceFilter
#' 
#' @seealso \code{\link{VarianceFilter}}
#'
#' @family getCandidate Methods
#'
#' @return A vector containing all genes with variances less than the target gene
#'
#' @export
#' 
#' @examples
#'
#' # Using the included Alzheimer's dataset, filter out only those transcription factors with variance
#' # within 50% of the variance of MEF2C
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' variance.filter <- VarianceFilter(mtx.assay = mtx.sub, targetGene = "MEF2C")
#' tfs <- getCandidates(variance.filter)

setMethod("getCandidates", "VarianceFilter",

          function(obj){

              # Retrive the extra arguments
              mtx <- obj@mtx.assay
              target.gene <- obj@targetGene
              var.size <- obj@varSize

              # Designate the target genes and tfs              
              tfs <- setdiff(rownames(mtx), target.gene)              
              tf.mtx <- mtx[-c(which(rownames(mtx) == target.gene)),]              
              target.mtx <- mtx[which(rownames(mtx) == target.gene),]
              
              # Find the variances              
              tf.var <- apply(tf.mtx,1,stats::var)              
              target.var <- stats::var(target.mtx)

              # Return only the genes with variances within the var.size of target gene variance
              var.idx <- which(tf.var > (1-var.size)*target.var & tf.var < (1+var.size)*target.var)
              tfs <- names(tf.var[var.idx])
              tf.vars <- tf.var[var.idx]

              return(list(tfs = tfs, tf.vars = tf.vars))
          }
)
#----------------------------------------------------------------------------------------------------
