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
#' @param regions A data frame that specifies the regions of interest
#' (default = data.frame())
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
#' genome.db.uri <- paste("sqlite:/",db.address,"mef2c.neighborhood.hg38.gtfAnnotation.db", sep = "/")
#' project.db.uri <- paste("sqlite:/",db.address,"mef2c.neigborhood.hg38.footprints.db", sep = "/")
#' target.gene <- "MEF2C"
#' size.upstream <- 1000
#' size.downstream <- 1000
#'
#' # Construct a Trena object and use it to retrive the regions
#' trena <- Trena("hg38")
#' regions <- getProximalPromoter(trena,target.gene, size.upstream, size.downstream)
#'  
#' footprint.filter <- FootprintFilter(genomeDB = genome.db.uri, footprintDB = project.db.uri,
#' regions = regions)

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
#' # in the included Alzheimer's dataset, using the Trena object to get regions
#' target.gene <- "MEF2C"
#' db.address <- system.file(package="trena", "extdata")
#' genome.db.uri <- paste("sqlite:/",db.address,"mef2c.neighborhood.hg38.gtfAnnotation.db", sep = "/")
#' project.db.uri <- paste("sqlite:/",db.address,"mef2c.neigborhood.hg38.footprints.db", sep = "/")
#' size.upstream <- 1000
#' size.downstream <- 1000
#'
#' # Construct a Trena object and use it to retrive the regions
#' trena <- Trena("hg38")
#' regions <- getProximalPromoter(trena,target.gene, size.upstream, size.downstream)
#'
#' footprint.filter <- FootprintFilter(genomeDB = genome.db.uri, footprintDB = project.db.uri,
#' regions = regions)
#' footprints <- getCandidates(footprint.filter)

setMethod("getCandidates", "FootprintFilter",

          function(obj){
              
              # Retrieve the FootprintFinder object and find the footprints
              fp <- FootprintFinder(obj@genomeDB, obj@footprintDB)
#              tbl.out <- data.frame()

              # Assume the regions are all different genes; put them into a list
              tbl.out <- list()
              for(r in 1:nrow(obj@regions)){
                  chromLoc <- as.list(obj@regions[r,])
                  if(!obj@quiet) printf(" FootprintFilter::getCandidates, getFootprintsInRegion %s-%s",
                                        chromLoc$start, chromLoc$end)
                  tbl.fp <- try(with(chromLoc, getFootprintsInRegion(fp, chrom, start, end)))
                  # If it's an error, return that error
                  if(class(tbl.fp) == "try-error"){
                      warning("FootprintFinder error with region %s-%s",
                              chromLoc$start, chromLoc$end)
                      tbl.fp <- "No footprints found"
                  } else {
                      # Rename the "name" column to "motifName"
                      names(tbl.fp)[which(names(tbl.fp) == "name")] <- "motifName"
                  }
                  #                  tbl.out <- rbind(tbl.out, tbl.fp)
                  tbl.out[[r]] <- tbl.fp
              } # for region
              
              # Close the DB connections and return the tabl
              closeDatabaseConnections(fp)
              return(tbl.out)
          }) # getCandidates
#----------------------------------------------------------------------------------------------------
