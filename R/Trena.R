#----------------------------------------------------------------------------------------------------
#' @name Trena-class
#' @rdname Trena-class
#' @aliases Trena
#'
#' @import methods

.Trena <- setClass ("Trena",
                    representation = representation(
                        genomeName="character",
                        quiet="logical")
                    )
#----------------------------------------------------------------------------------------------------
setGeneric('getRegulatoryChromosomalRegions',  signature='obj',
           function(obj, chromosome, chromStart, chromEnd, regulatoryRegionSources, targetGene, targetGeneTSS,
                    combine=FALSE) standardGeneric("getRegulatoryChromosomalRegions"))

setGeneric('getRegulatoryTableColumnNames',  signature='obj', function(obj) standardGeneric ('getRegulatoryTableColumnNames'))

setGeneric('getGeneModelTableColumnNames',  signature='obj', function(obj) standardGeneric ('getGeneModelTableColumnNames'))

setGeneric('createGeneModel', signature='obj', function(obj, targetGene,  solverNames, tbl.regulatoryRegions, mtx)
    standardGeneric('createGeneModel'))

setGeneric('getProximalPromoter', signature='obj', function(obj, geneSymbol, tssUpstream, tssDownstream)
    standardGeneric('getProximalPromoter'))
#----------------------------------------------------------------------------------------------------
# a temporary hack: some constants
genome.db.uri <- "postgres://bddsrds.globusgenomics.org/hg38"   # has gtf and motifsgenes tables
#----------------------------------------------------------------------------------------------------
#' Define an object of class Trena
#'
#' @description
#' The Trena class provides a convenient wrapper for the most commonly used filters and solver in the \code{trena}
#' package. Given a particular genome (one of \code{c("hg38","mm10")}, the Trena class provides methods to
#' retrieve information about possible regulators for a target gene, assess the effects of SNPs, and create gene models
#' using the flexible \code{\link{EnsembleSolver}} class.
#'
#' @rdname Trena-class
#'
#' @param genomeName A string indicating the genome used by the Trena object. Currently, only human and mouse (
#' "hg38","mm10") are supported
#' @param quiet A logical indicating whether or not the Trena object should print output
#'
#' @return An object of the Trena class
#'
#' @export
#'
#' @examples
#' # Create a Trena object using the human hg38 genome
#' trena <- Trena("hg38")
#'
#' @seealso \code{\link{getRegulatoryChromosomalRegions}}, \code{\link{getRegulatoryTableColumnNames}},
#' \code{\link{getGeneModelTableColumnNames}}, \code{\link{createGeneModel}}

Trena = function(genomeName, quiet=TRUE)
{
    stopifnot(genomeName %in% c("hg38", "mm10"))

    obj <- .Trena(genomeName=genomeName, quiet=quiet)

    obj

} # constructor
#----------------------------------------------------------------------------------------------------
#' Retrieve the column names in the regulatory table for a Trena object
#'
#' @rdname getRegulatoryTableColumnNames
#' @aliases getRegulatoryTableColumnNames
#'
#' @param obj An object of class Trena
#'
#' @return A character vector listing the column names in the Trena object regulatory table
#'
#' @export
#'
#' @examples
#' # Create a Trena object and retrieve the column names of the regulatory table
#' trena <- Trena("mm10")
#' tbl.cols <- getRegulatoryTableColumnNames(trena)

setMethod('getRegulatoryTableColumnNames', 'Trena',

          function(obj){
              c("chrom", "motifStart", "motifEnd", "motifName", "strand", "score", "length", "distance.from.tss", "id")
          })

#----------------------------------------------------------------------------------------------------
#' Retrieve the column names in the gene model table for a Trena object
#'
#' @rdname getGeneModelTableColumnNames
#' @aliases getGeneModelTableColumnNames
#'
#' @param obj An object of class Trena
#'
#' @return A character vector listing the column names in the Trena object gene model table
#'
#' @export
#'
#' @examples
#' # Create a Trena object and retrieve the column names of the gene model table
#' trena <- Trena("mm10")
#' tbl.cols <- getRegulatoryTableColumnNames(trena)

setMethod('getGeneModelTableColumnNames', 'Trena',

          function(obj){
              c("tf", "randomForest", "pearson", "spearman", "betaLasso", "pcaMax", "concordance")
          })

#----------------------------------------------------------------------------------------------------
.callFootprintFilter <- function(obj, source, chromosome, chromStart, chromEnd, targetGene, targetGeneTSS)
{
    chromLocString <- sprintf("%s:%d-%d", chromosome, chromStart, chromEnd)
    fpFilter <- FootprintFilter(genome.db.uri, source,
                                regions=data.frame(chrom=chromosome, start=chromStart, end=chromEnd, stringsAsFactors=FALSE))
    tbl.fp <- getCandidates(fpFilter)

    if(nrow(tbl.fp) == 0){
        warning("no footprints found in %s:%d-%d, targetGene is %s", chromosome, chromStart, chromEnd, targetGene);
        return(tbl.fp)
    }

    tbl.fp <- tbl.fp[, c("chrom", "start", "endpos", "motifName", "length", "strand", "score1", "score2", "score3")]
    colnames(tbl.fp) <- c("chrom", "motifStart", "motifEnd", "motifName", "length", "strand", "score1", "score", "score3")

    distance <- tbl.fp$motifStart - targetGeneTSS
    direction <- rep("upstream", length(distance))
    direction[which(distance < 0)] <- "downstream"
    tbl.fp$distance.from.tss <- distance
    tbl.fp$id <- sprintf("%s.fp.%s.%06d.%s", targetGene, direction, abs(distance), tbl.fp$motifName)
    # a final rearrangement of columns, to match our standard
    tbl.fp <- tbl.fp[, getRegulatoryTableColumnNames(obj)]

    tbl.fp

} # .callFootprintFilter
#----------------------------------------------------------------------------------------------------
.callHumanDHSFilter <- function(obj, chromosome, chromStart, chromEnd, targetGene, targetGeneTSS)
{
    if(!obj@quiet) printf("--- in .callHumanDHS")
    chromLocString <- sprintf("%s:%d-%d", chromosome, chromStart, chromEnd)
    dhsFilter <- HumanDHSFilter(genome="hg38",
                                encodeTableName="wgEncodeRegDnaseClustered",
                                pwmMatchPercentageThreshold=85L,
                                geneInfoDatabase.uri=genome.db.uri,
                                regions=data.frame(chrom=chromosome,
                                                   start=chromStart,
                                                   end=chromEnd,
                                                   stringsAsFactors=FALSE),
                                pfms = as.list(query(query(MotifDb, "sapiens"),
                                                     "jaspar2016"))
                                )

    tbl.dhs <- getCandidates(dhsFilter)
    if(nrow(tbl.dhs) == 0)
        return(tbl.dhs)

    tbl.dhs$length <- nchar(tbl.dhs$match)
    distance <- tbl.dhs$motifStart - targetGeneTSS
    direction <- rep("upstream", length(distance))
    direction[which(distance < 0)] <- "downstream"

    colnames(tbl.dhs)[grep("motifRelativeScore", colnames(tbl.dhs))] <- "score"
    #colnames(tbl.dhs)[grep("tfs", colnames(tbl.dhs))] <- "tf"
    tbl.dhs$distance.from.tss <- distance
    tbl.dhs$id <- sprintf("%s.dhs.%s.%06d.%s", targetGene, direction, abs(distance), tbl.dhs$motifName)

    tbl.dhs <- tbl.dhs[, getRegulatoryTableColumnNames(obj)]

    tbl.dhs

} # .callHumanDHSFilter
#----------------------------------------------------------------------------------------------------
#' Get the regulatory chromosomal regions for a Trena object
#'
#' @rdname getRegulatoryChromosomalRegions
#' @aliases getRegulatoryChromosomalRegions
#'
#' @param obj An object of class Trena
#' @param chromosome A choromosome of interest
#' @param chromStart The beginning of the desired region
#' @param chromEnd The end of the desired region
#' @param regulatoryRegionSources A vector containing the names of sources for chromosome information. These can be
#' addresses of footprint databases or the names of DHS databases
#' @param targetGene A target gene of interest
#' @param targetGeneTSS An integer giving the location of the target gene's transcription start site
#' @param combine A logical indicating whether or not to combine the output into one data frame (default = FALSE)
#'
#' @export
#'
#' @return A list of regulatory regions for the supplied target gene. If \code{combine} is set to \code{TRUE},
#' the list is converted into a data frame.
#'
#' @examples
#' # Get regulatory regions for MEF2C from a footprint database
#' database.filename <- system.file(package="trena", "extdata", "mef2c.neigborhood.hg38.footprints.db")
#' database.uri <- sprintf("sqlite://%s", database.filename)
#' sources <- c(database.uri)
#'
#' trena <- Trena("hg38")
#' chromosome <- "chr5"
#' mef2c.tss <- 88904257
#' loc.start <- mef2c.tss - 1000
#' loc.end   <- mef2c.tss + 1000
#'
#' regions <- getRegulatoryChromosomalRegions(trena, chromosome, mef2c.tss-1000, mef2c.tss+1000,
#' sources, "MEF2C", mef2c.tss)
#'
#' # Get regulatory regions for AQP4 from a Human DHS source
#' trena <- Trena("hg38")
#' aqp4.tss <- 26865884
#' chromosome <- "chr18"
#' sources <- c("encodeHumanDHS")
#'
#' regions <- getRegulatoryChromosomalRegions(trena, chromosome, aqp4.tss-1, aqp4.tss+3, sources, "AQP4", aqp4.tss)

setMethod('getRegulatoryChromosomalRegions', 'Trena',

          function(obj, chromosome, chromStart, chromEnd, regulatoryRegionSources, targetGene, targetGeneTSS,
                   combine=FALSE){

              tbl.combined <- data.frame()
              result <- list()
              # some bookeeeping to permit duplicate sources, useful only in testing
              source.count <- 0
              all.source.names <- regulatoryRegionSources

              encodeDHS.source.index <- grep("encodeHumanDHS", regulatoryRegionSources)

              if(length(encodeDHS.source.index)){
                  source.count <- source.count + 1
                  regulatoryRegionSources <- regulatoryRegionSources[-encodeDHS.source.index]
                  printf("calling HumanDHSFilter, span: %d",  1 + chromEnd - chromStart);
                  tbl.dhs <- .callHumanDHSFilter(obj, chromosome, chromStart, chromEnd, targetGene, targetGeneTSS)
                  result[[source.count]] <- tbl.dhs
                  if(combine)
                      tbl.combined <- rbind(tbl.combined, tbl.dhs)
              } # if encode DSH source requested


              for(source in regulatoryRegionSources){
                  source.count <- source.count + 1
                  printf("calling footprintFilter with source = '%s', span: %d", source, 1 + chromEnd - chromStart);
                  tbl.fp <- .callFootprintFilter(obj, source, chromosome, chromStart, chromEnd, targetGene, targetGeneTSS)
                  if(combine)
                      tbl.combined <- rbind(tbl.combined, tbl.fp)
                  result[[source.count]] <- tbl.fp
              } # for source
              names(result) <- all.source.names
              if(combine)
                  result[["all"]] <- tbl.combined
              result
          }) # getRegulatoryChromosomalRegions

#----------------------------------------------------------------------------------------------------
#' Create a model for a target gene using a Trena object
#'
#' @rdname createGeneModel
#' @aliases createGeneModel
#'
#' @param obj An object of class Trena
#' @param targetGene The name of a target gene to use for building a model
#' @param solverNames A character vector containing the solver names to be used for building the model
#' @param tbl.regulatoryRegions A data frame of regulatory regions, typically generated by using a filter
#' @param mtx An assay matrix of expression data
#'
#' @return A data frame containing the gene model
#'
#' @export
#'
#' @examples
#' # Create a Trena object for human and make a gene model for "MEF2C" using a footprint filter
#' trena <- Trena("hg38")
#' chromosome <- "chr5"
#' mef2c.tss <- 88904257
#' loc.start <- mef2c.tss - 1000
#' loc.end   <- mef2c.tss + 1000
#'
#' database.filename <- system.file(package="trena", "extdata", "mef2c.neigborhood.hg38.footprints.db")
#' database.uri <- sprintf("sqlite://%s", database.filename)
#' sources <- c(database.uri)
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#'
#' motifs.list <- getRegulatoryChromosomalRegions(trena, chromosome, mef2c.tss-1000, mef2c.tss+1000,
#' sources, "MEF2C", mef2c.tss)
#'
#' library(MotifDb)
#' tbl.motifs.tfs <- associateTranscriptionFactors(MotifDb, motifs.list[[1]], source="MotifDb", expand.rows=TRUE)
#' model.mef2c <- createGeneModel(trena, "MEF2C", c("lasso","ridge","randomforest"), tbl.motifs.tfs, mtx.sub)


setMethod('createGeneModel', 'Trena',

          function(obj, targetGene, solverNames, tbl.regulatoryRegions, mtx){

              stopifnot("geneSymbol" %in% colnames(tbl.regulatoryRegions))
              unique.tfs.from.regulatory.regions <- unique(tbl.regulatoryRegions$geneSymbol)
              tfs <- intersect(unique.tfs.from.regulatory.regions, rownames(mtx))
              if(!obj@quiet)
                  printf("tf candidate count, in mtx, in tbl.regulatory.regions: %d/%d", length(tfs),
                         length(unique.tfs.from.regulatory.regions))

              if(length(tfs) == 0)
                  return(data.frame())

              solver <- EnsembleSolver(mtx, targetGene=targetGene, candidateRegulators=tfs, solverNames)
              tbl.model <- run(solver)
              tbl.tf.frequencies <- as.data.frame(table(tbl.regulatoryRegions$geneSymbol))
              colnames(tbl.tf.frequencies) <- c("gene", "bindingSites")
              tbl.model <- merge(tbl.model, tbl.tf.frequencies, by="gene")
              tbl.model <- tbl.model[order(tbl.model$pcaMax, decreasing=TRUE),]
              tbl.model
          }) # createGeneModel
#----------------------------------------------------------------------------------------------------
#' Grab the region of the proximal promoter for a given gene symbol
#'
#' For the genome of a given Trena object, retrieve a data frame containing the region
#' surrounding a target gene.
#'
#' @rdname getProximalPromoter
#' @aliases getProximalPromoter
#'
#' @param obj An object of class Trena
#' @param geneSymbol A gene of interest
#' @param tssUpstream A designated distance upstream of the promoter to use as a shoulder
#' (default = 1000)
#' @param tssDownstream A designated distance downstream of the promoter to use as a shoulder
#' (default = 1000)
#'
#' @return A dataframe containing the regions surrounding the proximal promoter
#'
#' @export
#'
#' @examples
#' # Retrieve the proximal promoter for MEF2C using a shoulder size of 2000 on each side
#' trena <- Trena("hg38")
#' regions <- getProximalPromoter(trena, "MEF2C", 2000, 2000)

setMethod("getProximalPromoter", "Trena",

          function(obj, geneSymbol,
                   tssUpstream = 1000,
                   tssDownstream = 1000){

              # Switch the name of the database and filter we use
              db.name <- switch(obj@genomeName,
                                "hg38" = "hsapiens_gene_ensembl",
                                "mm10" = "mmusculus_gene_ensembl")
              filter.name <- switch(obj@genomeName,
                                "hg38" = "hgnc_symbol",
                                "mm10" = "mgi_symbol")

              my.mart <- biomaRt::useMart(biomart="ensembl", dataset= db.name)              
              
              tbl.geneInfo <- biomaRt::getBM(attributes=c("chromosome_name",
                                                          "transcription_start_site",
                                                          filter.name,
                                                          "strand"),
                                             filters=filter.name, value=geneSymbol, mart=my.mart)
              
              if(nrow(tbl.geneInfo) == 0)
                  return(NA)
              
              # make sure all transcripts are on the same strand
              strand <- unique(tbl.geneInfo$strand)
              stopifnot(length(strand) == 1)
              chrom <- sprintf("chr%s", unique(tbl.geneInfo$chromosome_name))
              stopifnot(length(chrom) == 1)
              
              # assume + strand
              tss <- min(tbl.geneInfo$transcription_start_site)
              start.loc <- tss - tssUpstream
              end.loc   <- tss + tssDownstream
              
              if(strand == -1){
                  tss <- max(tbl.geneInfo$transcription_start_site)
                  start.loc <- tss - tssDownstream
                  end.loc   <- tss + tssUpstream
              }
                    
          return (data.frame(chrom=chrom, start=start.loc, end=end.loc, stringsAsFactors=FALSE))

          })# getProximalPromoter
#----------------------------------------------------------------------------------------------------
