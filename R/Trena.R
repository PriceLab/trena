#' @import MotifDb
#' @importFrom DBI   dbConnect dbListTables dbGetQuery dbListConnections dbDisconnect
#' @importFrom RPostgreSQL dbConnect dbListTables dbGetQuery dbListConnections dbDisconnect
#' @importFrom stats fivenum prcomp sd
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

setGeneric('createGeneModelFromRegulatoryRegions', signature='obj',
           function(obj, targetGene,  solverNames, tbl.regulatoryRegions, mtx)
              standardGeneric('createGeneModelFromRegulatoryRegions'))

setGeneric('createGeneModelFromTfList', signature='obj',
           function(obj, targetGene,  solverNames, tfList, mtx)
              standardGeneric('createGeneModelFromTfList'))

setGeneric('getProximalPromoter', signature='obj', function(obj, geneSymbols, tssUpstream, tssDownstream)
    standardGeneric('getProximalPromoter'))

setGeneric('assessSnp', signature='obj', function(obj, pfms, variant, shoulder, pwmMatchMinimumAsPercentage, relaxedMatchDelta=25)
    standardGeneric('assessSnp'))
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
#' \code{\link{getGeneModelTableColumnNames}}

Trena = function(genomeName, quiet=TRUE)
{
    stopifnot(genomeName %in% c("hg38", "mm10", "tair10"))

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
              c("tf", "randomForest", "pearson", "spearman", "betaLasso")
          })

#----------------------------------------------------------------------------------------------------
.callFootprintFilter <- function(obj, source, chromosome, chromStart, chromEnd, targetGene, targetGeneTSS)
{
    chromLocString <- sprintf("%s:%d-%d", chromosome, chromStart, chromEnd)
    fpFilter <- FootprintFilter(genome.db.uri, source,
                                regions=data.frame(chrom=chromosome, start=chromStart, end=chromEnd, stringsAsFactors=FALSE))
    list.fp <- getCandidates(fpFilter)
    tbl.fp <- do.call("rbind",list.fp)

    if(nrow(tbl.fp) == 0){
       warning(sprintf("no footprints found in %s, region %s:%d-%d, targetGene is %s",
                       source, chromosome, chromStart, chromEnd, targetGene));
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
    pfms <- as.list(query(query(MotifDb, "sapiens"), "jaspar2018"))
    dhsFilter <- HumanDHSFilter(genomeName="hg38",
                                encodeTableName="wgEncodeRegDnaseClustered",
                                pwmMatchPercentageThreshold=85L,
                                geneInfoDatabase.uri=genome.db.uri,
                                regions=data.frame(chrom=chromosome,
                                                   start=chromStart,
                                                   end=chromEnd,
                                                   stringsAsFactors=FALSE),
                                pfms = pfms
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
#' @rdname createGeneModelFromRegulatoryRegions
#' @aliases createGeneModelFromRegulatoryRegions
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
#' if(interactive()){  # takes too long for the bioconductor build
#'    # Create a Trena object for human and make a gene model for "MEF2C" using a footprint filter
#'    trena <- Trena("hg38")
#'    chromosome <- "chr5"
#'    mef2c.tss <- 88904257
#'    loc.start <- mef2c.tss - 1000
#'    loc.end   <- mef2c.tss + 1000
#'
#'    database.filename <- system.file(package="trena", "extdata", "mef2c.neigborhood.hg38.footprints.db")
#'    database.uri <- sprintf("sqlite://%s", database.filename)
#'    sources <- c(database.uri)
#'    load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#'
#'    motifs.list <- getRegulatoryChromosomalRegions(trena, chromosome, mef2c.tss-1000, mef2c.tss+1000,
#'    sources, "MEF2C", mef2c.tss)
#'
#'    library(MotifDb)
#'    tbl.motifs.tfs <- associateTranscriptionFactors(MotifDb, motifs.list[[1]], source="MotifDb", expand.rows=TRUE)
#'    model.mef2c <- createGeneModelFromRegulatoryRegions(trena, "MEF2C", c("lasso","ridge","randomforest"),
#'                                                        tbl.motifs.tfs, mtx.sub)
#'    } # if interactive


setMethod('createGeneModelFromRegulatoryRegions', 'Trena',

          function(obj, targetGene, solverNames, tbl.regulatoryRegions, mtx){

              stopifnot(is.data.frame(tbl.regulatoryRegions))
              stopifnot("geneSymbol" %in% colnames(tbl.regulatoryRegions))
              unique.tfs.from.regulatory.regions <- unique(tbl.regulatoryRegions$geneSymbol)
              tfs <- intersect(unique.tfs.from.regulatory.regions, rownames(mtx))
              if(!obj@quiet)
                  printf("tf candidate count, in mtx, in tbl.regulatory.regions: %d/%d", length(tfs),
                         length(unique.tfs.from.regulatory.regions))

              if(length(tfs) == 0)
                  return(data.frame())

             solver <- EnsembleSolver(mtx, targetGene=targetGene, candidateRegulators=tfs, solverNames,
                                      geneCutoff=1)
              tbl.model <- run(solver)
              tbl.tf.frequencies <- as.data.frame(table(tbl.regulatoryRegions$geneSymbol))
              colnames(tbl.tf.frequencies) <- c("gene", "bindingSites")
              tbl.model <- merge(tbl.model, tbl.tf.frequencies, by="gene")
              if("pearsonCoeff" %in% colnames(tbl.model))  # our default ordering
                 tbl.model <- tbl.model[order(abs(tbl.model$pearsonCoeff), decreasing=TRUE),]
              tbl.model
          }) # createGeneModelFromRegulatoryRegions

#----------------------------------------------------------------------------------------------------
#' Create a model for a target gene using a Trena object
#'
#' @rdname createGeneModelFromTfList
#' @aliases createGeneModelFromTfList
#'
#' @param obj An object of class Trena
#' @param targetGene The name of a target gene to use for building a model
#' @param solverNames A character vector containing the solver names to be used for building the model
#' @param tfList A character list, often the gene symbols for known transcription factors
#' @param mtx An assay matrix of expression data
#'
#' @return A data frame containing the gene model
#'
#' @export
#'
#' @examples
#' if(interactive()){  # takes too long for the bioconductor build
#'    # Create a Trena object for human and make a gene model for "MEF2C" using a footprint filter
#'    trena <- Trena("hg38")
#'    chromosome <- "chr5"
#'    mef2c.tss <- 88904257
#'    loc.start <- mef2c.tss - 1000
#'    loc.end   <- mef2c.tss + 1000
#'
#'    database.filename <- system.file(package="trena", "extdata", "mef2c.neigborhood.hg38.footprints.db")
#'    database.uri <- sprintf("sqlite://%s", database.filename)
#'    sources <- c(database.uri)
#'    load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#'
#'    model.mef2c <- createGeneModelFromTfList(trena, "MEF2C", c("lasso","ridge","randomforest"),
#'                                             tfList=c())
#'    } # if interactive


setMethod('createGeneModelFromTfList', 'Trena',

     function(obj, targetGene, solverNames, tfList, mtx){

        tfs <- unique(intersect(tfList, rownames(mtx)))
        if(!obj@quiet)
           printf("tf candidate count, also in mtx, from tfList: %d/%d", length(tfs), length(tfList))
        if(length(tfs) == 0)
            return(data.frame())
        solver <- EnsembleSolver(mtx, targetGene=targetGene, candidateRegulators=tfs, solverNames,
                                    geneCutoff=1)
        tbl.model <- run(solver)
        tbl.model$bindingSites <- NA
        if("rfScore" %in% colnames(tbl.model)){
            tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
        } else if ("pearsonCoeff" %in% colnames(tbl.model)){
           tbl.model <- tbl.model[order(abs(tbl.model$pearsonCoeff), decreasing=TRUE),]
        }


        tbl.model
        }) # createGeneModelFromRegulatoryRegions

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
#' @param geneSymbols A vector containing genes of interest
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
#' if(interactive()) {  # too slow for the bioc windows build
#'    # Retrieve the proximal promoter for MEF2C using a shoulder size of 2000 on each side
#'    trena <- Trena("hg38")
#'    regions <- getProximalPromoter(trena, "MEF2C", 2000, 2000)
#'    }

setMethod("getProximalPromoter", "Trena",

          function(obj, geneSymbols,
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
                                                          "transcript_tsl",
                                                          filter.name),
                                             filters=filter.name, value=geneSymbols, mart=my.mart)

              if(nrow(tbl.geneInfo) == 0)
                  return(NA)

              # Sort by hgnc_symbol and transcript_tsl, then pull the first entry for each gene
              tbl.geneInfo <- tbl.geneInfo[order(tbl.geneInfo[[filter.name]],
                                                 tbl.geneInfo$transcript_tsl),]
              tbl.geneInfo <- tbl.geneInfo[match(unique(tbl.geneInfo[[filter.name]]),
                                                 tbl.geneInfo[[filter.name]]),]

              # remove contigs and check to make sure it's just 1 chromosome
              chromosome_name <- NULL
              tbl.geneInfo <- subset(tbl.geneInfo, chromosome_name %in% c(1:22, "X", "Y", "MT"))
              chrom <- sprintf("chr%s", tbl.geneInfo$chromosome_name)

              tss <- tbl.geneInfo$transcription_start_site
              start.loc <- tss - tssDownstream
              end.loc   <- tss + tssUpstream

              return (data.frame(geneSymbol=tbl.geneInfo[[filter.name]],
                                 chrom=chrom,
                                 start=start.loc,
                                 end=end.loc,
                                 stringsAsFactors=FALSE))

          })# getProximalPromoter
#----------------------------------------------------------------------------------------------------
#' Assess the effect of a SNP using a Trena object
#'
#' @rdname assessSnp
#' @aliases assessSnp
#'
#' @param obj An object of class Trena
#' @param pfms A set of motif matrices, generally retrieved using MotifDb
#' @param variant A variant of interest
#' @param shoulder A distance from the TSS to use as a window
#' @param pwmMatchMinimumAsPercentage A minimum match percentage for the  motifs
#' @param relaxedMatchDelta A numeric indicating the degree of the match (default = 25)
#'
#' @return A data frame containing the gene model
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a Trena object for human, assign a variant, then assess the effects of the variant
#' trena <- Trena("hg38")
#'
#' library(MotifDb)
#' jaspar.human.pfms <- as.list(query(query(MotifDb, "jaspar2016"), "sapiens"))[21:25]
#'
#' variant <- "rs3875089" # chr18:26865469  T->C
#'
#' tbl <- assessSnp(trena, jaspar.human.pfms, variant, shoulder = 3,
#' pwmMatchMinimumAsPercentage = 65)
#' }

setMethod('assessSnp', 'Trena',

          function(obj, pfms, variant, shoulder, pwmMatchMinimumAsPercentage, relaxedMatchDelta=25){

              motifMatcher <- MotifMatcher(genomeName=obj@genomeName, pfms=pfms, quiet=obj@quiet)
              tbl.variant <- try(.parseVariantString(motifMatcher, variantString=variant), silent=TRUE)
              if(is(tbl.variant, "try-error")){
                  printf("error, unrecognized variant name: '%s'", variant)
                  return(data.frame())
              }
              tbl.regions <- data.frame(chrom=tbl.variant$chrom,
                                        start=tbl.variant$loc-shoulder,
                                        end=tbl.variant$loc+shoulder,
                                        stringsAsFactors=FALSE)

              tbl.wt  <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions,
                                                        pwmMatchMinimumAsPercentage=pwmMatchMinimumAsPercentage)
              if(nrow(tbl.wt) == 0){
                  warning(sprintf("no motifs found in reference sequence in neighborhood of %s with shoulder %d",
                                  variant, shoulder))
                  return(data.frame())
              }

              tbl.mut <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions,
                                                        pwmMatchMinimumAsPercentage=pwmMatchMinimumAsPercentage,
                                                        variants=variant)
              if(nrow(tbl.mut) == 0){
                  warning(sprintf("no motifs altered by %s with shoulder %d", variant, shoulder))
                  return(data.frame())
              }

              tbl.wt$signature <- sprintf("%s;%s;%s", tbl.wt$motifName, tbl.wt$motifStart, tbl.wt$strand)
              tbl.mut$signature <- sprintf("%s;%s;%s", tbl.mut$motifName, tbl.mut$motifStart, tbl.mut$strand)

              # comine wt and mut tables, reorder columns and rows for easier comprehension
              tbl <- rbind(tbl.wt[, c(1,12,2,3,4,5,6,7,8,13, 14)], tbl.mut[, c(1,12,2,3,4,5,6,7,8,13, 14)])
              tbl <- tbl[order(tbl$motifName, tbl$motifRelativeScore, decreasing=TRUE),]
              #tbl$signature <- sprintf("%s;%s;%s", tbl$motifName, tbl$motifStart, tbl$strand)
              #tbl <- tbl[,c(1,2,3:110)]

              # now look for less stringent matches.  these will be matched up with the
              # wt and mut motifs which do not yet have partners, thus enabling us to
              # provide a wt->mut motifScore.delta for each
              relaxedMatchPercentage <- pwmMatchMinimumAsPercentage-relaxedMatchDelta
              tbl.wt.relaxed <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, relaxedMatchPercentage)
              tbl.wt.relaxed$signature <- sprintf("%s;%s;%s", tbl.wt.relaxed$motifName, tbl.wt.relaxed$motifStart, tbl.wt.relaxed$strand)
              tbl.mut.relaxed <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, relaxedMatchPercentage, variants=variant)
              tbl.mut.relaxed$signature <- sprintf("%s;%s;%s", tbl.mut.relaxed$motifName, tbl.mut.relaxed$motifStart, tbl.mut.relaxed$strand)

              status <- NULL   # workaround "no global binding" error
              signatures.in.both <- intersect(subset(tbl, status=="mut")$signature, subset(tbl, status=="wt")$signature)
              signatures.only.in.wt <- setdiff(subset(tbl, status=="wt")$signature, subset(tbl, status=="mut")$signature)
              signatures.only.in.mut <- setdiff(subset(tbl, status=="mut")$signature, subset(tbl, status=="wt")$signature)

              tbl$assessed <- rep("failed", nrow(tbl))

              if(length(signatures.in.both) > 0) {
                 indices <- sort(unlist(lapply(signatures.in.both, function(sig) grep(sig, tbl$signature))))
                 tbl$assessed[indices] <- "in.both"
                 }

              if(length(signatures.only.in.wt) > 0) {
                 indices <- sort(unlist(lapply(signatures.only.in.wt, function(sig) grep(sig, tbl$signature))))
                 tbl$assessed[indices] <- "wt.only"
                 }

              if(length(signatures.only.in.mut) > 0) {
                 indices <- sort(unlist(lapply(signatures.only.in.mut, function(sig) grep(sig, tbl$signature))))
                 tbl$assessed[indices] <- "mut.only"
                 }

              tbl$delta <- 0

                 # find the mut scores for each of the "wt.only" entries, subtract from the wt score
                 # workaround "no visible binding"
              assessed <- NULL
              motifRelativeScore <- NULL

              tbl.wt.only  <- subset(tbl, assessed=="wt.only", select=c(signature, motifRelativeScore))
              if(nrow(tbl.wt.only) > 0){
                  sigs <- tbl.wt.only$signature
                  tbl.mut.scores <- subset(tbl.mut.relaxed, signature %in% sigs, select=c(signature, motifRelativeScore))
                  deltas <- unlist(lapply(sigs, function(sig){wt.score  <- subset(tbl.wt.only, signature==sig)$motifRelativeScore;
                      mut.score <- subset(tbl.mut.scores, signature==sig)$motifRelativeScore;
                      delta <- wt.score - mut.score
                  }))
                  tbl$delta[match(sigs, tbl$signature)] <- deltas
              } # if some wt.only entries

              # find the wt scores for each of the "mut.only" entries, subtract from the mut score

              tbl.mut.only  <- subset(tbl, assessed=="mut.only", select=c(signature, motifRelativeScore))
              if(nrow(tbl.mut.only) > 0){
                  sigs <- tbl.mut.only$signature
                  # find the wt scores for these muts, looking in the relaxedMatchPercentage match table
                  tbl.wt.scores <- subset(tbl.wt.relaxed, signature %in% sigs, select=c(signature, motifRelativeScore))
                  deltas <- unlist(lapply(sigs, function(sig){mut.score  <- subset(tbl.mut.only, signature==sig)$motifRelativeScore;
                      wt.score <- subset(tbl.wt.scores, signature==sig)$motifRelativeScore;
                      delta <- wt.score - mut.score
                  }))
                  tbl$delta[match(sigs, tbl$signature)] <- deltas
              } # tbl.mut.only > 0

              coi <-  c("motifName", "status", "assessed", "motifRelativeScore", "delta",
                        "signature", "chrom", "motifStart", "motifEnd", "strand", "match")

              stopifnot(all(coi %in% colnames(tbl)))
              tbl <- tbl[, coi]
              tbl$variant <- variant
              tbl
          }) # assessSnp
#----------------------------------------------------------------------------------------------------
