#' @title Create a FootprintFilter object
#'
#' @description
#' A FootprintFilter object allows for filtering based on gene footprinting databases. Using its
#' associated \code{getCandidates} method and URIs for both a genome database and project database,
#' a FootprintFilter object can be used to filter a list of possible transcription factors to those
#' that match footprint motifs for a supplied target gene.
#'
#' @include CandidateFilter.R
#' @include FootprintFinder.R
#' @import methods
#'
#' @name FootprintFilter-class
#' @rdname FootprintFilter-class
#' @aliases FootprintFilter

#----------------------------------------------------------------------------------------------------
.FootprintFilter <- setClass("FootprintFilter", contains = "CandidateFilter",
                             slots=list(genomeDB="character",
                                        footprintDB="character",
                                        regions="data.frame",
                                        footprintFinder="FootprintFinder")
                             )

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
#' @rdname FootprintFilter-class
#'
#' @param genomeDB A connection to a database that contains genome information
#' @param footprintDB A connection to a database that contains footprint information
#' @param geneCenteredSpec A named list that specifies the regions of interest using a target gene
#' and upstream/downstream sizes. If specified,it t should contain the following fields
#' (default = list()):
#' \itemize{
#' \item{"targetGene" A designated target gene that should be part of the mtx.assay data}
#' \item{"tssUpstream" An integer denoting the distance upstream of the target gene
#' to look for footprints}
#' \item{"tssDownstream" An integer denoting the distance downstream of the target gene
#' to look for footprints}
#' }
#' @param regionsSpec A character vector that specifies the regions of interest directly, using a
#' string containing chromosome number, starting position, and ending position. The string should
#' be formatted as follows: "chr##:start-end" (e.g. "chr1:10000-20000"). (default = list())
#' @param quiet A logical denoting whether or not the filter should print output
#'
#' @seealso \code{\link{getCandidates-FootprintFilter}}
#'
#' @export
#'
#' @return An object of the FootprintFilter class
#'
#' @family Filtering Objects
#'
#' @examples
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' db.address <- system.file(package="trena", "extdata")
#' genome.db.uri <- paste("sqlite:/",db.address,"genome.sub.db", sep = "/")
#' project.db.uri <- paste("sqlite:/",db.address,"project.sub.db", sep = "/")
#' target.gene <- "MEF2C"
#' size.upstream <- 1000
#' size.downstream <- 1000
#' geneCenteredSpec <- list(targetGene = target.gene, tssUpstream = size.upstream,
#' tssDownstream = size.downstream)
#' footprint.filter <- FootprintFilter(genomeDB = genome.db.uri, footprintDB = project.db.uri,
#' geneCenteredSpec = geneCenteredSpec)

FootprintFilter <- function(genomeDB, footprintDB, regions=data.frame(), quiet=TRUE)
{

   .FootprintFilter(CandidateFilter(quiet = quiet),
                    genomeDB=genomeDB,
                    footprintDB=footprintDB,
                    regions=regions)

} # FootprintFilter, the constructor
#----------------------------------------------------------------------------------------------------
#' Get candidate genes using the footprint filter
#'
#' @aliases getCandidates-FootprintFilter
#'
#' @param obj An object of class FootprintFilter
#'
#' @seealso \code{\link{FootprintFilter}}
#'
#' @family getCandidate Methods
#'
#' @return A list, where one element is the transcription factors found in the footprints and the
#' other is a data frame containing all the meta data for the footprints
#'
#' @export
#'
#' @examples
#'
#' # Use footprint filter with the included SQLite database for MEF2C to filter candidates
#' # in the included Alzheimer's dataset with a gene-centered spec
#' target.gene <- "MEF2C"
#' db.address <- system.file(package="trena", "extdata")
#' genome.db.uri <- paste("sqlite:/",db.address,"genome.sub.db", sep = "/")
#' project.db.uri <- paste("sqlite:/",db.address,"project.sub.db", sep = "/")
#' size.upstream <- 1000
#' size.downstream <- 1000
#' geneCenteredSpec <- list(targetGene = target.gene, tssUpstream = size.upstream, tssDownstream = size.downstream)
#' footprint.filter <- FootprintFilter(genomeDB = genome.db.uri, footprintDB = project.db.uri,
#' geneCenteredSpec = geneCenteredSpec)
#' tfs <- getCandidates(footprint.filter)
#'
#' # Perform the same operation, but use a region spec
#' mef2c.tss <- 88904257 # Empirically assign the MEF2C TSS as the location
#' chrom <- "chr5"
#' start <- mef2c.tss - 1000
#' end <- mef2c.tss + 1000
#' regionsSpec <- sprintf("%s:%d-%d", chrom, start, end)
#' footprint.filter <- FootprintFilter(genomeDB = genome.db.uri,
#' footprintDB = project.db.uri, regionsSpec = regionsSpec)
#' tfs <- getCandidates(footprint.filter)


setMethod("getCandidates", "FootprintFilter",

     function(obj){

         # Retrieve the FootprintFinder object and find the footprints
         fp <- FootprintFinder(obj@genomeDB, obj@footprintDB)
         tbl.out <- data.frame()

        #for(region in obj@regions){
           #chromLoc <- parseChromLocString(region)
        for(r in 1:nrow(obj@regions)){
           chromLoc <- as.list(obj@regions[r,])
           if(!obj@quiet) printf(" FootprintFilter::getCandidates, getFootprintsInRegion %s", region)
           tbl.fp <- try(with(chromLoc, getFootprintsInRegion(fp, chrom, start, end)))
           if(class(tbl.fp) == "try-error"){
              warning("FootprintFinder error with region %s", region)
              closeDatabaseConnections(fp)
              return(NULL)
              }
           tbl.out <- rbind(tbl.out, tbl.fp)
           } # for region

        closeDatabaseConnections(fp)
        tbl.out
        }) # getCandidates

#----------------------------------------------------------------------------------------------------
