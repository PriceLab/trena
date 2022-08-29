#' @title FeatureTable
#' @description A template for building documented, tested R6 classes
#' @name FeatureTable

#'
#' @examples
#'   rt <- FeatureTable()
#'
#' import R6
#' import data.table
#'
#' @export

FeatureTable = R6Class("FeatureTable",

    #--------------------------------------------------------------------------------
    private = list(tbl=NULL,
                   reference.genome=NULL,
                   target.gene=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(
         #' @description
         #' Creates a new instance of this [R6][R6::R6Class] class.
         #' @param target.gene character, a recognized gene symbol
         #' @param reference.genome character, a recognized genome shorthand code, e.g., "hg38"
         #' @return a new instance o FeatureTable
        initialize = function(target.gene, reference.genome){
            private$target.gene <- target.gene
            private$reference.genome <- reference.genome
            private$tbl <- as.data.table(data.frame())
            },
        #------------------------------------------------------------
        #' @description accessor for the table
        #' @return the data.table
        getTable = function(){
            private$tbl
            },

        #------------------------------------------------------------
        #' @description establish the basic genomic regions, typically
        #'   transcription factor binding sites
        #' @param tbl.regions a data.frame with at least chrom, start end columns
        setFundamentalRegions = function(tbl.regions){
            stopifnot(all(c("chrom", "start", "end") %in% colnames(tbl.regions)))
            private$tbl <- as.data.table(tbl.regions)
            },

        #------------------------------------------------------------
        #' @description annotate each of the intersecting fundamental regions with
        #'   the value of this feature
        #' @param tbl.bed a data.frame with chrom, start and end columns, and values
        #'        of interest described in the annotation guidel
        #' @param feature.genome character, a recognized genome shorthand code, e.g., "hg38"

        #' @param annotation.guide a neme list identifying columns of interest in
        #'       tbl.bed, and the column names to use for them in the feature table
        addFeature = function(tbl.feature, feature.genome, feature.guide){
            stopifnot(all(c("chrom", "start", "end") %in% colnames(tbl.feature)))
            stopifnot(feature.genome == self$reference.genome)
            tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.feature), GRanges(private$tbl)))
            if(nrow(tbl.ov) == 0) return()
            for(feature in names(feature.guide)){
                source.feature.name <- feature.guide[[feature]]
                feature.class <- class(tbl.feature[[source.feature.name]])
                vec <- vector(feature.class, length=nrow(private$tbl))
                vec[tbl.ov$subjectHits] <- tbl.feature[[source.feature.name]][tbl.ov$queryHits]
                private$tbl[[feature]] <- vec
                } # for feature
            } # addFeature

       ) # public

    ) # class
#--------------------------------------------------------------------------------
