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
#' @name HumanDHSFilter-class
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
                                       geneCenteredSpec="list",
                                       regionsSpec="character",
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
setGeneric("getSequence_tmp", signature="obj", function(obj, tbl.regions) standardGeneric ("getSequence_tmp"))
setGeneric("geneSymbolToTSS", signature="obj", function(obj, geneSymbol) standardGeneric("geneSymbolToTSS"))
#----------------------------------------------------------------------------------------------------
#' Create a CandidateFilter using Human DNAse Hypersensitivity
#'
#' @param genomeName A character string indicating the reference genome; currently, the only
#' accepted strings are "hg38" and "hg19", both of which are human genomes. 
#' @param encodeTableName (default = "wgEncodeRegDnaseClustered")
#' @param pwmMatchPercentageThreshold A numeric from 0-100 to serve as a threshold for a match
#' @param geneInfoDatabase.uri An address for a gene database
#' @param geneCenteredSpec A gene-centered spec
#' @param regionsSpec A regions spec
#' @param variants Variants
#' @param quiet A logical denoting whether or not the solver should print output
#'
#' @return A Solver class object with Random Forest as the solver
#'
#' @seealso  \code{\link{solve.RandomForest}}, \code{\link{getAssayData}}
#'
#' @family Solver class objects
#'
#' @export
#'
#' @examples
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' targetGene <- "MEF2C"
#' candidateRegulators <- setdiff(rownames(mtx.sub), targetGene)
#' rf.solver <- RandomForestSolver(mtx.sub, targetGene, candidateRegulators)

HumanDHSFilter <- function(genomeName,
                           encodeTableName="wgEncodeRegDnaseClustered",
                           #fimoDatabase.uri,
                           pwmMatchPercentageThreshold,
                           geneInfoDatabase.uri,
                           geneCenteredSpec=list(),
                           regionsSpec=list(),
                           variants=NA_character_,
                           quiet=TRUE)
{
    regions <- c();   # one or more chromLoc strings: "chr5:88903257-88905257"

    uri <- "http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt"
    x <- .readRawJasparMatrices(uri)
         # normalize, so that a frequency sum of 1.0 is true across the 4 possible bases at each position
    pfms <- lapply(x, function(e) apply(e$matrix, 2, function(col) col/sum(col)))
    names(pfms) <- as.character(lapply(x, function(e) e$title))

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
                   #fimoDB=fimo.db,
                   pwmMatchPercentageThreshold=pwmMatchPercentageThreshold,
                   encodeTableName=encodeTableName,
                   geneInfoDatabase.uri=geneInfoDatabase.uri,
                   genome=reference.genome,
                   geneCenteredSpec=geneCenteredSpec,
                   regionsSpec=regionsSpec,
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
#' #

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

setMethod("show", "HumanDHSFilter",

     function(object){
        s <- sprintf("HumanDHSFilter...")
        cat(s, sep="\n")
        })

#----------------------------------------------------------------------------------------------------
#' Get candidate genes using a human DHS filter
#'
#' @aliases getCandidates-FootprintFilter
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
#' goFilter <- GeneOntologyFilter(org.Hs.eg.db, GOTerm="GO:0006351")
#' candidates <- getCandidates(goFilter)

setMethod("geneSymbolToTSS", "HumanDHSFilter",

     function(obj, geneSymbol){
        geneInfo.db.info <- .parseDatabaseUri(obj@geneInfoDatabase.uri)
        host <- geneInfo.db.info$host
        dbname <- geneInfo.db.info$name
        driver <- RPostgreSQL::PostgreSQL()
        db.geneInfo <- DBI::dbConnect(driver, user= "trena", password="trena", dbname=dbname, host=host)
        #print(dbListTables(db.geneInfo))
        #db.gtf <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="gtf", host="whovian")
        query <- sprintf("select * from hg38human where moleculetype='gene' and gene_biotype='protein_coding' and gene_name='%s'",
                         obj@geneCenteredSpec$targetGene)
        tbl <- dbGetQuery(db.geneInfo, query);
        tss <- tbl$start[1];
        chrom <- tbl$chr
        DBI::dbDisconnect(db.geneInfo)
        list(chrom=chrom, tss=tss)
        })

#----------------------------------------------------------------------------------------------------
#' Get candidate genes using a human DHS filter
#'
#' @aliases getCandidates-FootprintFilter
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
#' goFilter <- GeneOntologyFilter(org.Hs.eg.db, GOTerm="GO:0006351")
#' candidates <- getCandidates(goFilter)

setMethod("getCandidates", "HumanDHSFilter",

    function(obj){

       regions <- c();   # one or more chromLoc strings: "chr5:88903257-88905257"

         # if present, extract a region from the geneCenteredSpec
       if(length(obj@geneCenteredSpec) == 3){
          expected.fields <- c("targetGene", "tssUpstream", "tssDownstream")
          stopifnot(all(expected.fields %in% names(obj@geneCenteredSpec)))
          x <- geneSymbolToTSS(obj)
          start <- x$tss - obj@geneCenteredSpec$tssUpstream
          end <- x$tss + obj@geneCenteredSpec$tssDownstream
          new.region.chromLocString <- sprintf("%s:%d-%d", x$chrom, start, end)
          regions <- c(regions, new.region.chromLocString)
          }

          # are there any explicit regions to include?
       if(!(all(is.na(obj@regionsSpec)))){
          regions <- c(regions, obj@regionsSpec)
          }

       printf("HumanDHSFilter::getCandidates");
       tbl.regions <- as.data.frame(do.call("rbind", lapply(regions, function(region) .parseChromLocString(region))))
       tbl.regions$chrom <- as.character(tbl.regions$chrom)
       tbl.regions$start <- as.integer(tbl.regions$start)
       tbl.regions$end <- as.integer(tbl.regions$end)
       tbl.dhs <- data.frame()
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
       mm <- MotifMatcher(genomeName=obj@genomeName, pfms=jaspar.human.pfms)
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
#' Get candidate genes using a human DHS filter
#'
#' @aliases getCandidates-FootprintFilter
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
#' goFilter <- GeneOntologyFilter(org.Hs.eg.db, GOTerm="GO:0006351")
#' candidates <- getCandidates(goFilter)

setMethod("getRegulatoryRegions", "HumanDHSFilter",

    function(obj, encode.table.name, chromosome, start, end, score.threshold=0) {

       driver <- RMySQL::MySQL()
       host <- "genome-mysql.cse.ucsc.edu"
       user <- "genome"
       dbname <- obj@genomeName

       db <- dbConnect(driver, user = user, host = host, dbname = dbname)

       #schema <- colnames(dbGetQuery(db, sprintf("select * from %s limit 1", tableName)))
       #suppressWarnings(dbGetQuery(db, sprintf("select * from %s limit 3", tableName)))

       #if(!all(c("chrom", "chromStart", "chromEnd") %in% schema)){
       #   printf("%s lacks chrom, start, end", tableName)
       #   lapply(dbListConnections(driver), dbDisconnect)
       #   return(data.frame())
       #   }

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

   #tbl.regions$motif.start <- -1 + tbl.regions$regionStart + motif.start
   #tbl.regions$motif.end <- -1 + tbl.regions$regionEnd + motif.end
   invisible(tbl.regions[, c("chrom", "chromStart", "chromEnd",  "count", "score")])

   }) # getRegulatoryRegions

#----------------------------------------------------------------------------------------------------
#' Get candidate genes using a human DHS filter
#'
#' @aliases getCandidates-FootprintFilter
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
#' goFilter <- GeneOntologyFilter(org.Hs.eg.db, GOTerm="GO:0006351")
#' candidates <- getCandidates(goFilter)

setMethod("getSequence_tmp", "HumanDHSFilter",

   function(obj, tbl.regions){
     gr.regions <- with(tbl.regions, GRanges(seqnames=chrom, IRanges(start=chromStart, end=chromEnd)))
     seqs <- BSgenome::getSeq(obj@genome, gr.regions)
     as.character(seqs)
     })  # getSequence_tmp

#----------------------------------------------------------------------------------------------------
