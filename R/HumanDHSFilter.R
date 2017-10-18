#' @title Create a HumanDHSFilter object
#'
#' @description
#' A HumanDHSFilter object allows for filtering based on DNAse hypersensitivity (DHS) data. Its
#' associated \code{getCandidates} method uses a genome from a BSgenome database (either hg19 or
#' hg38), DNA region specifications, and (variants/pfms,encodetablename, match%)
#' to filter a list of possible regulators factors to those that match the supplied criteria.
#'
#' @include CandidateFilter.R
#' @import methods
#' @import BSgenome
#'
#' @rdname HumanDHSFilter-class
#' @aliases HumanDHSFilter

#----------------------------------------------------------------------------------------------------
.HumanDHSFilter <- setClass("HumanDHSFilter",
                            contains="CandidateFilter",
                            slots=list(genomeName="character",
                                       genome="BSgenome",
                                       encodeTableName="character",
                                       pwmMatchPercentageThreshold="numeric",
                                       geneInfoDatabase.uri="character",   # access to gtf database
                                       regions="data.frame",
                                       variants="character",
                                       pfms="list",
                                       quiet="logical"))
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
setGeneric("getEncodeRegulatoryTableNames", signature="obj", function(obj) standardGeneric ("getEncodeRegulatoryTableNames"))

setGeneric("getRegulatoryRegions", signature="obj",
           function(obj, encode.table.name, chromosome, start, end, score.threshold=200, quiet=TRUE)
               standardGeneric ("getRegulatoryRegions"))
#----------------------------------------------------------------------------------------------------
#' Create a CandidateFilter using Human DNAse Hypersensitivity
#'
#' @param genomeName A character string indicating the reference genome; currently, the only
#' accepted strings are "hg38" and "hg19", both of which are human genomes.
#' @param encodeTableName (default = "wgEncodeRegDnaseClustered")
#' @param pwmMatchPercentageThreshold A numeric from 0-100 to serve as a threshold for a match
#' @param geneInfoDatabase.uri An address for a gene database
#' @param regions A data frame containing the regions of interest
#' @param variants A character vector containing a list of variants
#' @param pfms A list of position frequency matrices, often converted from a MotifList object
#' created by a MotifDb query
#' @param quiet A logical denoting whether or not the solver should print output
#'
#' @return A CandidateFilter class object that filters using Human DHS data
#'
#' @seealso  \code{\link{getCandidates-HumanDHSFilter}},
#'
#' @family Solver class objects
#'
#' @export
#'
#' @rdname HumanDHSFilter-class
#'
#' @examples
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' targetGene <- "VRK2"
#' promoter.length <- 1000
#' genomeName <- "hg38"
#' db.address <- system.file(package="trena", "extdata")
#' genome.db.uri    <- paste("sqlite:/", db.address, "vrk2.neighborhood.hg38.gtfAnnotation.db",  sep = "/")
#' 
#' # Grab regions for VRK2 using shoulder size of 1000
#' trena <- Trena(genomeName)
#' tbl.regions <- getProximalPromoter(trena, "VRK2", 1000, 1000)
#'
#' hd.filter <- HumanDHSFilter(genomeName, pwmMatchPercentageThreshold = 85,
#' geneInfoDatabase.uri = genome.db.uri, regions = tbl.regions,
#' pfms = as.list(query(query(MotifDb, "sapiens"),"jaspar2016")))

HumanDHSFilter <- function(genomeName,
                           encodeTableName="wgEncodeRegDnaseClustered",
                           pwmMatchPercentageThreshold,
                           geneInfoDatabase.uri,
                           regions,                           
                           variants=NA_character_,
                           pfms,
                           quiet=TRUE)
{   
    if(genomeName == "hg38"){
        # library(BSgenome.Hsapiens.UCSC.hg38) ## Remove the library reference
        reference.genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    }
    else if(genomeName == "hg19"){
        # library(BSgenome.Hsapiens.UCSC.hg19) ## Remove the library reference
        reference.genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    }
    else {
        stop(sprintf("HumanDHSFilter genome.name not in hg19, hg38: '%s'", genomeName))
    }
    
    .HumanDHSFilter(CandidateFilter(quiet = quiet),
                    genomeName=genomeName,
                    pwmMatchPercentageThreshold=pwmMatchPercentageThreshold,
                    encodeTableName=encodeTableName,
                    geneInfoDatabase.uri=geneInfoDatabase.uri,
                    genome=reference.genome,
                    regions=regions,
                    variants=variants,
                    pfms=pfms,
                    quiet=quiet)
    
} # HumanDHSFilter, the constructor
#----------------------------------------------------------------------------------------------------
#' Get Encode regulatory tables using a human DHS filter
#'
#' @rdname getEncodeRegulatoryTableNames-HumanDHSFilter
#' @aliases getEncodeRegulatoryTableNames
#'
#' @param obj An object of class HumanDHSFilter
#'
#' @seealso \code{\link{HumanDHSFilter}}
#'
#' @return A character vector containing the names of the Encode regulatory tables for the regions
#' contained in the HumanDHSFilter object
#'
#' @export
#'
#' @examples
#'
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' targetGene <- "VRK2"
#' promoter.length <- 1000
#' genomeName <- "hg38"
#' db.address <- system.file(package="trena", "extdata")
#' genome.db.uri    <- paste("sqlite:/", db.address, "vrk2.neighborhood.hg38.gtfAnnotation.db",  sep = "/")
#' jaspar.human <- as.list(query(query(MotifDb, "sapiens"),"jaspar2016"))
#'
#' # Grab regions for VRK2 using shoulder size of 1000
#' trena <- Trena(genomeName)
#' tbl.regions <- getProximalPromoter(trena, "VRK2", 1000, 1000)
#'
#' hd.filter <- HumanDHSFilter(genomeName, pwmMatchPercentageThreshold = 85,
#' geneInfoDatabase.uri = genome.db.uri, regions = tbl.regions, pfms = jaspar.human)
#'
#' getEncodeRegulatoryTableNames(hd.filter) 

setMethod("getEncodeRegulatoryTableNames", "HumanDHSFilter",

          function(obj){
              driver <- RMySQL::MySQL()
              host <- "genome-mysql.cse.ucsc.edu"
              user <- "genome"
              dbname <- obj@genomeName
              db <- DBI::dbConnect(driver, user = user, host = host, dbname = dbname)
              all.tableNames <- DBI::dbListTables(db);
              # manual check (13 apr 2017) shows that only wgEncodeReg Peak tabel
              # and the "wgEncodeRegDnaseClustered" table, have scored chromosomal regions in them
              tableNames <- grep("wgEncodeReg.*Peak$", DBI::dbListTables(db), value=TRUE)
              clusteredTable <- switch(obj@genomeName,
                                       hg19="wgEncodeRegDnaseClusteredV3",
                                       hg38="wgEncodeRegDnaseClustered",
                                       NA)
              if(clusteredTable %in% all.tableNames)
                  tableNames <- c(clusteredTable, tableNames)
              lapply(dbListConnections(driver), DBI::dbDisconnect)
              invisible(tableNames)
          })
#----------------------------------------------------------------------------------------------------
#' Show the details of a human DHS filter
#'
#' @aliases show-HumanDHSFilter
#'
#' @param object An object of class HumanDHSFilter
#'
#' @seealso \code{\link{HumanDHSFilter}}
#'
#' @return A list, where one element a character vector of transcription factors that match
#' the GO term and the other is an empty data frame.
#'
#' @export
#'
#' @examples
#' # Make a filter and show it
#' #' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' targetGene <- "VRK2"
#' promoter.length <- 1000
#' genomeName <- "hg38"
#' db.address <- system.file(package="trena", "extdata")
#' genome.db.uri    <- paste("sqlite:/", db.address, "vrk2.neighborhood.hg38.gtfAnnotation.db",  sep = "/")
#' jaspar.human <- as.list(query(query(MotifDb, "sapiens"),"jaspar2016"))
#'
#' # Grab regions for VRK2 using shoulder size of 1000
#' trena <- Trena(genomeName)
#' tbl.regions <- getProximalPromoter(trena, "VRK2", 1000, 1000)
#'
#' hd.filter <- HumanDHSFilter(genomeName, pwmMatchPercentageThreshold = 85,
#' geneInfoDatabase.uri = genome.db.uri, regions = tbl.regions, pfms = jaspar.human)
#'
#' show(hd.filter)

setMethod("show", "HumanDHSFilter",

          function(object){
              s <- sprintf("HumanDHSFilter...")
              cat(s, sep="\n")
          })
#----------------------------------------------------------------------------------------------------
#' Get candidate genes using a human DHS filter
#'
#' @aliases getCandidates-HumanDHSFilter
#'
#' @param obj An object of class FootprintFilter
#'
#' @seealso \code{\link{FootprintFilter}}
#'
#' @family getCandidate Methods
#'
#' @return A list, where one element a character vector of transcription factors that match
#' the GO term and the other is an empty data frame.
#'
#' @export
#'
#' @examples
#'
#' # Make a filter for "transcription, DNA-templated" and use it to filter candidates
#' #' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' targetGene <- "VRK2"
#' promoter.length <- 1000
#' genomeName <- "hg38"
#' db.address <- system.file(package="trena", "extdata")
#' genome.db.uri    <- paste("sqlite:/", db.address, "vrk2.neighborhood.hg38.gtfAnnotation.db",  sep = "/")
#' jaspar.human <- as.list(query(query(MotifDb, "sapiens"),"jaspar2016"))
#'
#' # Grab regions for VRK2 using shoulder size of 1000
#' trena <- Trena(genomeName)
#' tbl.regions <- getProximalPromoter(trena, "VRK2", 1000, 1000)
#'
#' hd.filter <- HumanDHSFilter(genomeName, pwmMatchPercentageThreshold = 85,
#' geneInfoDatabase.uri = genome.db.uri, regions = tbl.regions, pfms = jaspar.human)
#'
#' getCandidates(hd.filter)

setMethod("getCandidates", "HumanDHSFilter",

          function(obj){
              
              tbl.regions <- obj@regions

              if(!obj@quiet){
                  printf("HumanDHSFilter::getCandidates, from these regions:");
                  print(tbl.regions)
              }
              
              tbl.dhs <- data.frame()
              
              if(!obj@quiet){
                  printf("HumanDHSFilter::getCandidates, getRegulatoryRegions for %d regions", nrow(tbl.regions))
                  print(tbl.regions)
              }
              
              for(r in 1:nrow(tbl.regions)){
                  tbl.new <- getRegulatoryRegions(obj, obj@encodeTableName, tbl.regions$chrom[r], tbl.regions$start[r], tbl.regions$end[r])
                  tbl.dhs <- rbind(tbl.dhs, tbl.new)
              }
              if(!obj@quiet)
                  printf("found %d DHS regions in %d requested regions", nrow(tbl.dhs), nrow(tbl.regions))
              
              if(nrow(tbl.dhs) == 0){
                  return(tbl.dhs)
              }
              
              colnames(tbl.dhs) <- c("chrom", "start", "end", "count", "score")
              jaspar.human.pfms <- as.list(query(query(MotifDb, "hsapiens"), "jaspar2016"))
              mm <- MotifMatcher(genomeName=obj@genomeName, pfms=jaspar.human.pfms, obj@quiet)
              tbl <- findMatchesByChromosomalRegion(mm, tbl.dhs[, 1:3],
                                                    pwmMatchMinimumAsPercentage=obj@pwmMatchPercentageThreshold,
                                                    variants=obj@variants)
              if(!obj@quiet)
                  printf(" and %d motifs", nrow(tbl))
              
              preferred.colnames <- c("motifName", "chrom", "motifStart", "motifEnd", "strand", "motifScore", "motifRelativeScore",
                                      "match", "regulatoryRegionStart", "regualtoryRegionEnd", "regulatorySequence", "variant")
              colnames(tbl) <- preferred.colnames
              long.sequences <- which(tbl$regulatorySequence > 20)
              if(length(long.sequences) > 0){
                  tbl$regulatorySequence[long.sequences] <- paste(substr(tbl$regulatorySequence[long.sequences], 1, 17), "...", sep="")
              }
              tbl
          }) # getCandidates
#----------------------------------------------------------------------------------------------------
#' Get a tabel of regulatory regions for a Human DHS filter
#'
#' @aliases getRegulatoryRegions
#' @rdname getRegulatoryRegions
#'
#' @param obj An object of class HumanDHSFilter
#' @param encode.table.name A vector of names for Encode tables
#' @param chromosome The chromosome of interest
#' @param start The starting position
#' @param end The ending position
#' @param score.threshold A threshold for the score (default = 200)
#'
#' @seealso \code{\link{HumanDHSFilter}}
#'
#' @return A data frame containing the regulatory regions for the filter, including the
#' chromosome, start, and end positions, plus the count and score of each region.
#'
#' @export
#'
#' @examples
#'
#' # Make a filter for "transcription, DNA-templated" and use it to filter candidates
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' targetGene <- "VRK2"
#' promoter.length <- 1000
#' genomeName <- "hg38"
#' db.address <- system.file(package="trena", "extdata")
#' genome.db.uri    <- paste("sqlite:/", db.address, "vrk2.neighborhood.hg38.gtfAnnotation.db",  sep = "/")
#' jaspar.human <- as.list(query(query(MotifDb, "sapiens"),"jaspar2016"))
#'
#' # Grab regions for VRK2 using shoulder size of 1000
#' trena <- Trena(genomeName)
#' tbl.regions <- getProximalPromoter(trena, "VRK2", 1000, 1000)
#'
#' hd.filter <- HumanDHSFilter(genomeName, pwmMatchPercentageThreshold = 85,
#' geneInfoDatabase.uri = genome.db.uri, regions = tbl.regions, pfms = jaspar.human)
#'
#' chrom <- "chr2"
#' rs13384219.loc <- 57907323
#' start <- rs13384219.loc - 10
#' end <- rs13384219.loc + 10
#'
#' tableNames <- getEncodeRegulatoryTableNames(hd.filter)
#'
#' getRegulatoryRegions(hd.filter, tableNames[1], chrom, start, end)

setMethod("getRegulatoryRegions", "HumanDHSFilter",
          
          function(obj, encode.table.name, chromosome, start, end, score.threshold=0) {
              
              driver <- RMySQL::MySQL()
              host <- "genome-mysql.cse.ucsc.edu"
              user <- "genome"
              dbname <- obj@genomeName
              
              if(!obj@quiet)
                  printf("connecting to %s/%s/%s as %s", host, dbname, encode.table.name,  user);
              
              db <- dbConnect(driver, user = user, host = host, dbname = dbname)
              
              main.clause <- sprintf("select * from %s where", encode.table.name);
              
              # Pull out the regions corresponding to the region in ENCODE
              query <- paste(main.clause,
                             sprintf("chrom = '%s'", chromosome),
                             sprintf("and chromStart >= %d", start),
                             sprintf("and chromEnd <= %d", end),
                             collapse = " ")
              
              if(!obj@quiet)
                  printf("query: %s", query)
              
              # handle the usual case first: a start:end region many times larger than a typical DHS region
              suppressWarnings(  # MySQL returns unsigned integers.  hide these unproblematic conversion warnings
                  tbl.regions <- dbGetQuery(db, query)
              )
              
              if(!obj@quiet)
                  printf("%d DHS regions reported in %d bases, start:end unmodified", nrow(tbl.regions), 1 + end - start)
              
              # if no hits, then perhaps a very small region is requested, one which falls entirely within a DHS region
              if(nrow(tbl.regions) == 0) {
                  extension <- 10000
                  if(!obj@quiet)
                      printf("possible that start:end (%d) is small relative to DHS regions, extend by %d", (1 + end - start),
                             extension);
                  query <- paste(main.clause,
                                 sprintf("chrom = '%s'", chromosome),
                                 sprintf("and chromStart >= %d", start - extension),
                                 sprintf("and chromEnd   <= %d", end + extension),
                                 collapse = " ")
                  if(!obj@quiet)
                      printf("query with extension: %s", query)
                  suppressWarnings(  # MySQL returns unsigned integers.  hide these unproblematic conversion warnings
                      tbl.regionsExtended <- dbGetQuery(db, query)
                  )
                  if(!obj@quiet) printf("query with extended region, %d rows", nrow(tbl.regionsExtended))
                  if(nrow(tbl.regionsExtended) > 0) { # now find just the intersection of DHS and requested region
                      if(!obj@quiet)
                          printf("tbl.regionsExtended: %d rows, now have intersection", nrow(tbl.regionsExtended));
                      gr.regions <- with(tbl.regionsExtended,
                                         GRanges(seqnames=chromosome,
                                                 IRanges::IRanges(chromStart, chromEnd)))
                      gr.target <- GRanges(seqnames=chromosome, IRanges::IRanges(start, end))
                      # use range to collapse any multiple intersections down to that of the original target
                      gr.intersect <- GenomicRanges::intersect(gr.target, gr.regions, ignore.strand=TRUE)
                      if(!obj@quiet) printf("GenomicRanges intersections of extended region with original target: %d", length(gr.intersect))
                      if(length(gr.intersect) >= 1){
                          tbl.ov <- as.data.frame(findOverlaps(gr.intersect, gr.regions, type="any"))
                          tbl.regions <- cbind(as.data.frame(gr.intersect[tbl.ov$queryHits]),
                                               tbl.regionsExtended[tbl.ov$subjectHits, c("name", "score")])
                          colnames(tbl.regions) <- c("chrom", "chromStart", "chromEnd", "width", "strand", "name", "score")
                      } # one or more region from the extended query intersects with the requested start:end.
                  } # small region query, within DHS region
              } # small region, extension search
              
              lapply(dbListConnections(driver), dbDisconnect)
              
              tbl.regions$chrom <- as.character(tbl.regions$chrom)
              
              # the ucsc database call the  itemCount columns "name".  fix that
              if("name" %in% colnames(tbl.regions)){
                  colnames(tbl.regions)[match("name", colnames(tbl.regions))] <- "count"
              }
              
              invisible(tbl.regions[, c("chrom", "chromStart", "chromEnd",  "count", "score")])
              
          }) # getRegulatoryRegions
#----------------------------------------------------------------------------------------------------
