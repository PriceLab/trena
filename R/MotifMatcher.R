#' @title Create a MotifMatcher object
#'
#' @description
#' A MotifMatcher object is used directly by the \code{\link{HumanDHSFilter}} class to match motif
#' matrices to where they occur in the supplied genome.
#'
#' @import methods
#' @import BSgenome
#'
#' @rdname MotifMatcher-class
#' @aliases MotifMatcher
#----------------------------------------------------------------------------------------------------
.MotifMatcher <- setClass('MotifMatcher',
                          representation(genome="BSgenome",
                                         pfms="list",
                                         quiet="logical")
                          )
#----------------------------------------------------------------------------------------------------
setGeneric("getPfms",            signature="obj",
           function(obj) standardGeneric ("getPfms"))

setGeneric("getSequence",        signature="obj",
           function(obj,tbl.regions,variants=NA_character_)
               standardGeneric ("getSequence"))

setGeneric(".parseVariantString", signature="obj",
           function(obj, variantString)
               standardGeneric (".parseVariantString"))

setGeneric("findMatchesByChromosomalRegion", signature="obj",
           function(obj, tbl.regions, pwmMatchMinimumAsPercentage, variants=NA_character_)
               standardGeneric ("findMatchesByChromosomalRegion"))
#----------------------------------------------------------------------------------------------------
#' @title Class MotifMatcher
#' @name MotifMatcher-class
#' @rdname MotifMatcher-class
#'
#' @description
#' The MotifMatcher class is used to match motif position weight matrices to places where they occur
#' in a given genome. It requries specification of a genome to search in and a list of motifs to
#' search for. Ordinarily this class is primarily used by the HumanDHSFilter, but can alternatively
#' be used to search for motifs in a given genome without any filtering functionality.
#'
#' @param genomeName A character string identifying an object of type BSgenome. The genome object
#' contains the information for a specific human genome and should be either "hg38" or "hg19". The
#' supplied genome serves as the search space for matching motifs (default = "hg38").
#' @param pfms A list of motif matrices to serve as queries for the target genome. If supplied,
#' these should be created using a MotifList object from the \code{\link{MotifDb}} package (see
#' example below). If unspecified, the motifs will default to all vertebrates in the JASPAR
#' database (default = list())
#' @param quiet A logical denoting whether or not the MotifMatcher object should print output
#'
#' @return An object of the MotifMatcher class
#'
#' @export
#'
#' @seealso \code{\link{HumanDHSFilter}}
#'
#' @examples
#' # Specify the genome, and motif list to create a MotifMatcher for only human motifs
#' library(MotifDb)
#' mm <- MotifMatcher( genomeName="hg38",
#' pfms = as.list(query(MotifDb, "sapiens")))

MotifMatcher <- function(genomeName,
                         pfms,
                         quiet=TRUE)
{
    genomeName <- tolower(genomeName)

    stopifnot(is.list(pfms))
    stopifnot(genomeName %in% c("hg19", "hg38", "saccer3", "tair10"))

    if(genomeName == "hg38"){
       reference.genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
       }

    else if(genomeName == "hg19"){
       reference.genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
        }

    else if(genomeName == "mm10"){
       reference.genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
       }

    else if(tolower(genomeName) == "saccer3"){
       reference.genome <- BSgenome.Scerevisiae.UCSC.sacCer3::BSgenome.Scerevisiae.UCSC.sacCer3
       }

    else if(tolower(genomeName) == "tair10"){
       reference.genome <- BSgenome.Athaliana.TAIR.TAIR9::BSgenome.Athaliana.TAIR.TAIR9
       }

    else {
       stop(sprintf("MotifMatch, genomeName not in hg19, hg38: '%s'", genomeName))
       }

    .MotifMatcher(genome=reference.genome, pfms=pfms, quiet=quiet)

} # MotifMatcher constructor
#----------------------------------------------------------------------------------------------------
#' Show a MotifMatcher object
#'
#' @rdname show.MotifMatcher
#' @aliases show.MotifMatcher
#'
#' @param object An object of the class MotifMatcher
#'
#' @return A truncated view of the supplied object
#'
#' @examples
#' load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' lassopv.solver <- LassoPVSolver(mtx.sub, target.gene, tfs)
#' show(lassopv.solver)

setMethod("show", "MotifMatcher",

          function(object){
              s <- sprintf("MotifMatcher for genome %s", object@genome)
              cat(s, sep="\n")
          })
#----------------------------------------------------------------------------------------------------
#' Find Motif Matches by Chromosomal Region
#'
#' Given a MotifMatcher object, a table of chromosomal regions, and a minimum match percentage,
#' pull out a list containing a data frame of the motifs in those regions and a character vector
#' of their associated transcription factors.
#'
#' @rdname findMatchesByChromosomalRegion
#' @aliases findMatchesByChromosomalRegion
#'
#' @param obj An object of class MotifMatcher
#' @param tbl.regions A data frame where each row contains a chromosomal region with the fields
#' "chrom", "start", and "end".
#' @param pwmMatchMinimumAsPercentage A percentage (0-100) used as a cutoff for what constitutes
#' a motif match
#' @param variants A character containing variants to use for the matching (default = NA_character_).
#' The variants should either have the same number of entries as rows in the \code{tbl.regions},
#' or they should not be supplied.
#'
#' @return A list containing a data frame of the motifs in the given regions and a character
#' vector of their associated transcription factors
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform a simple match in the rs13384219 neighborhood
#' library(MotifDb)
#' motifMatcher <- MotifMatcher(genomeName="hg38",
#' pfms = as.list(query(query(MotifDb, "sapiens"),"jaspar2016")), quiet=TRUE)
#' tbl.regions <- data.frame(chrom="chr2", start=57907313, end=57907333, stringsAsFactors=FALSE)
#' x <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=92)
#'
#' # Perform the same match, but now include a variant
#' x.mut <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions,
#' pwmMatchMinimumAsPercentage=92, variants = "rs13384219")
#' }

setMethod("findMatchesByChromosomalRegion", "MotifMatcher",

          function(obj, tbl.regions, pwmMatchMinimumAsPercentage, variants=NA_character_){

              x <- lapply(1:nrow(tbl.regions),
                          function(r) getSequence(obj, tbl.regions[r,], variants))
              tbl.regions <- do.call(rbind, x)
              if(!obj@quiet){
                  printf("---- MotifMatcher::findMatchesByChromosomalRegion")
                  print(tbl.regions)
              }
              tbl.motifs.list <- .getScoredMotifs(tbl.regions$seq, obj@pfms,
                                                  pwmMatchMinimumAsPercentage, obj@quiet)

              region.count <- nrow(tbl.regions)
              tbl.out <- data.frame()
              all.tfs <- c()
              for(i in seq_len(region.count)){
                  tbl.motifs <- tbl.motifs.list[[i]]
                  if(nrow(tbl.motifs) > 0){
                      colnames(tbl.motifs) <- c("motifStart", "motifEnd", "width",
                                                "motifScore", "maxScore",
                                                "motifRelativeScore", "motifName",
                                                "match", "strand")
                      tbl.region <- tbl.regions[i,]
                      colnames(tbl.region) <- c("chrom", "chromStart", "chromEnd", "seq", "status")
                      tbl.out <- rbind(tbl.out, cbind(tbl.motifs, tbl.region))
                  }
              } # for i

              if(nrow(tbl.out) == 0)
                  return(tbl=data.frame())

              tbl.out$motifStart <- tbl.out$motifStart + tbl.out$chromStart - 1;
              tbl.out$motifEnd <- tbl.out$motifEnd  + tbl.out$chromStart - 1;
              # change some column names
              #colnames(tbl.out)[4] <- "motifscore"
              #colnames(tbl.out)[2] <- "endpos"
              #colnames(tbl.out)[7] <- "motifname"
              tbl.out <- tbl.out[, -c(3,5)] # get rid of width and maxScore columns

              desired.column.order <- c("motifName", "chrom", "motifStart", "motifEnd", "strand",
                                        "motifScore", "motifRelativeScore", "match",
                                        "chromStart", "chromEnd", "seq", "status") #, "count", "score")
              tbl.out <- tbl.out[, desired.column.order]
              tbl.out <- tbl.out[order(tbl.out$motifScore, decreasing=TRUE),]

                 # for MotifDb motif names, simplify, keeping only the final token:
                 # "Hsapiens-jaspar2016-FOXH1-MA0479.1" -> "MA0479.1"
              tokens <- strsplit(tbl.out$motifName, "-")
              short.motif.names <- unlist(lapply(tokens, function(tokenSet) tokenSet[length(tokenSet)]))
              tbl.out$shortMotif <- short.motif.names
                 # tbl.mg will soon come from MotifDb
                 #tbl.mg <- tfGeneSymbolForMotif(obj, tbl.out$motifName)
                 #tfs.by.motif <- lapply(tbl.out$motifName, function(m) subset(tbl.mg, motif==m)$geneSymbol)
                 #all.tfs <- sort(unique(unlist(tfs.by.motif)))
                 #tfs.by.motif.joined <- unlist(lapply(tfs.by.motif, function(m) paste(m, collapse=";")))
                 #tbl.out$tf <- tfs.by.motif.joined
              max.sequence.chars <- 20
              if(nchar(tbl.out$seq[1]) > max.sequence.chars)
                  tbl.out$seq <- paste(substring(tbl.out$seq, 1, max.sequence.chars-3), "...", sep="")
              return(tbl.out)

          }) #findMatchesByChromosomalRegion
#----------------------------------------------------------------------------------------------------
#' Retrieve the motifs from the pfms slot
#'
#' Given a MotifMatcher object, return the motifs, which are stored in the pfms slot.
#'
#' @rdname getPfms
#' @aliases getPfms
#'
#' @param obj An object of class MotifMatcher
#'
#' @return The list of motif matrices stored in the pfms slot.
#'
#' @export
#'
#' @examples
#'
#' # Return the default matrix of JASPAR motifs
#' library(MotifDb)
#' motifMatcher <- MotifMatcher(genomeName="hg38", pfms = as.list(query(MotifDb, "sapiens")))
#' motifs <- getPfms(motifMatcher)

setMethod("getPfms", "MotifMatcher",

          function(obj){
              return(obj@pfms)
          })
#----------------------------------------------------------------------------------------------------
.matchPwmForwardAndReverse <- function(sequence, pfm, motifName, min.match.percentage=95, quiet=TRUE)
{
   xyz <- .matchPwmForwardAndReverse

    min.match.as.string <- sprintf("%02d%%", min.match.percentage)

    hits.fwd <- Biostrings::matchPWM(pfm, sequence, with.score=TRUE, min.score=min.match.as.string)
    hits.rev <- Biostrings::matchPWM(Biostrings::reverseComplement(pfm),
                                     sequence, with.score=TRUE, min.score=min.match.as.string)

    max.score <-  Biostrings::maxScore(pfm)
    tbl <- data.frame()
    if(length(hits.fwd) > 0){
        if(!quiet) printf("%d +", length(hits.fwd))
        match <- substring(as.character(subject(hits.fwd)), start(ranges(hits.fwd)), end(ranges(hits.fwd)))
        actual.score <- mcols(hits.fwd)$score
        relative.score <- actual.score/max.score
        tbl <- data.frame(ranges(hits.fwd),
                          score=mcols(hits.fwd)$score, maxScore=max.score,relativeScore=relative.score,
                          motif=motifName, match=match, strand="+",
                          stringsAsFactors=FALSE)
    } # hits.fwd

    if(length(hits.rev) > 0){
        if(!quiet) printf("%d -", length(hits.rev))
        match <- substring(as.character(subject(hits.rev)), start(ranges(hits.rev)), end(ranges(hits.rev)))
        actual.score <- mcols(hits.rev)$score
        relative.score <- actual.score/max.score
        tbl.rev <- data.frame(ranges(hits.rev),
                              score=mcols(hits.rev)$score, maxScore=max.score, relativeScore=relative.score,
                              motif=motifName, match=match, strand="-",
                              stringsAsFactors=FALSE)
        # transform the start/end so that they are forward-strand relative
        #true.start <- 1 + nchar(sequence) - tbl.rev$end
        #true.end   <- 1 + nchar(sequence) - tbl.rev$start
        #tbl.rev$start <- true.start
        #tbl.rev$end   <- true.end
        tbl <- rbind(tbl, tbl.rev)
    }

    #printf("returning match fwd/bwd tbl with %d rows", nrow(tbl))
    #print(table(tbl$strand))

    tbl

} # .matchPwmForwardAndReverse
#----------------------------------------------------------------------------------------------------
.findMotifs <- function(sequence, pfms, min.match.percentage=95, quiet=TRUE)
{
    min.match.as.string <- sprintf("%02d%%", min.match.percentage)

    count <- length(pfms)
    xx <- lapply(1:count, function(i) {
        .matchPwmForwardAndReverse(sequence, pfms[[i]], names(pfms)[i], min.match.percentage, quiet)
    })

    tbl.result <- do.call("rbind", xx)
    if(nrow(tbl.result) == 0){
        return(data.frame())
    }
    else{
        tbl.result$motif <- as.character(tbl.result$motif)
        tbl.result$match <- as.character(tbl.result$match)
        tbl.result$strand <- as.character(tbl.result$strand)
        return(tbl.result[order(tbl.result$score,decreasing=TRUE),])
    }

}  # .findMotifs
#------------------------------------------------------------------------------------------------------------------------
.getScoredMotifs <- function(seqs, pfms, min.match.percentage=95, quiet=TRUE)
{
    parseLine <- function(textOfLine) {
        # first delete the leading A, C, G or T.  then the square brackets.  then convert
        x <- substr(textOfLine, 2, nchar(textOfLine))
        x2 <- sub(" *\\[ *", "", x)
        x3 <- sub(" *\\] *", "", x2)
        counts <- as.integer(strsplit(x3, "\\s+", perl=TRUE)[[1]])
        return(counts)
    } # parseLine

    parseJasparPwm = function (lines) {
        stopifnot(length(lines)==5) # title line, one line for each base
        motif.name.raw = strsplit(lines[1], "\t")[[1]][1]
        motif.name <- gsub(">", "", motif.name.raw, fixed=TRUE)
        # expect 4 rows, and a number of columns we can discern from  the incoming text.
        a.counts <- parseLine(lines[[2]])
        c.counts <- parseLine(lines[[3]])
        g.counts <- parseLine(lines[[4]])
        t.counts <- parseLine(lines[[5]])
        stopifnot(length(a.counts) == length(c.counts))
        stopifnot(length(a.counts) == length(g.counts))
        stopifnot(length(a.counts) == length(t.counts))
        cols <- length(a.counts)
        mtx <- matrix (nrow=4, ncol=cols, dimnames=list(c('A','C','G','T'), as.character(1:cols)))
        mtx[1,] <- a.counts
        mtx[2,] <- c.counts
        mtx[3,] <- g.counts
        mtx[4,] <- t.counts
        return(list(title=motif.name, matrix=mtx))
    } # parsePwm

    readRawJasparMatrices = function (uri) {
        all.lines <- scan(uri, what=character(0), sep='\n', quiet=TRUE)
        title.lines <- grep ('^>', all.lines)
        title.line.count <- length (title.lines)
        max <- title.line.count - 1
        pwms <- list()
        for(i in 1:max){
            start.line <- title.lines [i]
            end.line <- title.lines [i+1] - 1
            x <- parseJasparPwm (all.lines [start.line:end.line])
            pwms[[i]] <- list(title=x$title, matrix=x$matrix)
        } # for i
        invisible (pwms)
    } # readRawJasparMatrices

    xyz <- "MotifMatcher .getScoredMotifs"
    result <- lapply(seqs, function(seq) .findMotifs(seq, pfms, min.match.percentage, quiet))
    if(is.null(result))
        result <- data.frame()

    result

} # .getScoredMotifs
#----------------------------------------------------------------------------------------------------
.parseLine <- function(textOfLine)
{
    # first delete the leading A, C, G or T.  then the square brackets.  then convert
    x <- substr(textOfLine, 2, nchar(textOfLine))
    x2 <- sub(" *\\[ *", "", x)
    x3 <- sub(" *\\] *", "", x2)
    counts <- as.integer(strsplit(x3, "\\s+", perl=TRUE)[[1]])

    return(counts)

} # .parseLine
#------------------------------------------------------------------------------------------------------------------------
.parseJasparPwm = function (lines)
{
    stopifnot(length(lines)==5) # title line, one line for each base
    motif.name.raw = strsplit(lines[1], "\t")[[1]][1]
    motif.name <- gsub(">", "", motif.name.raw, fixed=TRUE)
    # expect 4 rows, and a number of columns we can discern from  the incoming text.
    a.counts <- .parseLine(lines[[2]])
    c.counts <- .parseLine(lines[[3]])
    g.counts <- .parseLine(lines[[4]])
    t.counts <- .parseLine(lines[[5]])
    stopifnot(length(a.counts) == length(c.counts))
    stopifnot(length(a.counts) == length(g.counts))
    stopifnot(length(a.counts) == length(t.counts))
    cols <- length(a.counts)
    mtx <- matrix (nrow=4, ncol=cols, dimnames=list(c('A','C','G','T'), as.character(1:cols)))
    mtx[1,] <- a.counts
    mtx[2,] <- c.counts
    mtx[3,] <- g.counts
    mtx[4,] <- t.counts

    return(list(title=motif.name, matrix=mtx))

} # .parsePwm
#------------------------------------------------------------------------------------------------------------------------
.readRawJasparMatrices = function (uri)
{
    all.lines <- scan(uri, what=character(0), sep='\n', quiet=TRUE)
    title.lines <- grep ('^>', all.lines)
    title.line.count <- length (title.lines)
    max <- title.line.count - 1
    pwms <- list()
    for(i in 1:max){
        start.line <- title.lines [i]
        end.line <- title.lines [i+1] - 1
        x <- .parseJasparPwm (all.lines [start.line:end.line])
        pwms[[i]] <- list(title=x$title, matrix=x$matrix)
    } # for i

    invisible (pwms)

} # .readRawJasparMatrices
#----------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------
# return a new version of tbl.regions:
#
#  1) with a "status" column added
#  2) in which regions with a variant have their sequence altered to include the alternate base.
#  3) if there # are multpiole alternate alleles, then a separate row is created for each alternate;
#  4) in which any regions without mutants  are returned as is (with the alt column left empty)
#
.injectSnp <- function(tbl.regions, tbl.variants)
{
    stopifnot(nrow(tbl.regions) == 1)
    stopifnot(all(c("chrom", "start", "end", "seq") %in% colnames(tbl.regions)))
    gr.regions  <- with(tbl.regions,  GRanges(seqnames=chrom, IRanges(start=start, end=end)))
    gr.variants <- with(tbl.variants, GRanges(seqnames=chrom, IRanges(start=loc,   end=loc)))
    suppressWarnings(
        tbl.overlaps <- as.data.frame(findOverlaps(gr.variants, gr.regions))
    )
    colnames(tbl.overlaps) <- c("variant", "region")

    if(nrow(tbl.overlaps) == 0)
        return(tbl.regions)

    regions.without.snps <- setdiff(1:nrow(tbl.regions), tbl.overlaps$region)
    tbl.regions.out <- tbl.regions[regions.without.snps,]
    tbl.regions.out$status <- rep("wt", nrow(tbl.regions.out))

    new.sequence <- tbl.regions$seq
    status.string <- ""

    for(r in seq_len(nrow(tbl.overlaps))){
        wt  <- as.list(tbl.regions[tbl.overlaps$region[r],])
        alt <- as.list(tbl.variants[tbl.overlaps$variant[r],])
        offset <- 1 + alt$loc - wt$start
        mut.loc <- wt$start + offset
        wt.base <- substring(new.sequence, offset, offset)
        new.sequence <- sprintf("%s%s%s", substr(new.sequence, 1, (offset-1)), alt$mut,
                                substr(new.sequence, offset+1, nchar(new.sequence)))
        #alt.string <- sprintf("%s:%d(%s->%s)", alt$chrom, alt$loc, alt$wt, alt$mut)
    }

    tbl.regions.out <- tbl.regions
    tbl.regions.out$seq <- new.sequence
    tbl.regions.out$status <- "mut"
    tbl.regions.out

    #tbl.regions.out[order(tbl.regions.out$start, decreasing=FALSE),]

} # .injectSnp
#----------------------------------------------------------------------------------------------------
# going obsolete.   see function .parseVariantString
setMethod(".parseVariantString", "MotifMatcher",

          function(obj, variantString) {
              # quick sanity check to detect obvious errors in the variant string
              is.rsid <- grepl("^rs", variantString)
              is.explicit.chromLoc <- grepl("^chr", variantString) &&
                  length(strsplit(variantString, ":")[[1]]) == 4
              is.plausible.variant.string <- is.rsid || is.explicit.chromLoc
              if(!is.plausible.variant.string)
                  stop(sprintf("variant '%s' is neither a plausible rsID nor an explicit chrom:loc:ref:var", variantString))
              if(grepl("^rs", variantString)){
                  explicitly.specified.alternate.allele <- NA
                  tokens <- strsplit(variantString, ":")[[1]]
                  if(length(tokens) == 2){  # eg, "rs3763040:A" when the snp includes multiple alternate alleles
                      variantString <- tokens[1]
                      explicitly.specified.alternate.allele <- tokens[2]
                  }
                  snp.info <- as.data.frame(BSgenome::snpsById(
                      SNPlocs.Hsapiens.dbSNP150.GRCh38::SNPlocs.Hsapiens.dbSNP150.GRCh38, variantString))[1,]
                  chrom <- as.character(snp.info$seqnames)
                  if(!grepl("ch", chrom))
                      chrom <- sprintf("chr%s", chrom)
                  if(!grepl("chr", chrom))
                      chrom <- sub("ch", "chr", chrom)
                  start <- as.numeric(snp.info$pos)
                  ambiguity.code <- snp.info$alleles_as_ambig
                  elements.string <- Biostrings::IUPAC_CODE_MAP[[ambiguity.code]]
                  elements <- strsplit(elements.string,'')[[1]]
                  wt <- as.character(BSgenome::getSeq(obj@genome, chrom, start, start))
                  mut <- setdiff(elements, wt)
                  snp <- list(chrom=chrom, loc=start, wt=wt, mut=mut)
                  variant.count <- length(snp$mut)
                  result <- vector(mode="list", length=variant.count)
                  for(v in 1:variant.count){
                      result[[v]] <- snp
                      result[[v]]$mut <- snp$mut[v]
                  } # for v
                  tbl.out <- do.call(rbind.data.frame, result)
                  if(variant.count > 1){
                      chosen.row <- 1  # the default if alternate allele not specified explicitly
                      if(is.na(explicitly.specified.alternate.allele)){
                          warning(sprintf("alternate allele of %d options not specified, choosing first: %s",
                                          variant.count, tbl.out[chosen.row, "mut"]))
                      }
                      else{
                          chosen.row <- match(explicitly.specified.alternate.allele, tbl.out$mut)
                          if(is.na(chosen.row)) stop(sprintf("alternate allele %s not in %s", tokens[2], variantString))
                      }
                      tbl.out <- tbl.out[chosen.row,]
                  } # variant.count > 1
              } # if rsid
              else{  # no rsid string supplied: instead, an explicit chrom:loc:wt:mut string ("chr2:57907323:A:G")
                  tokens <- strsplit(variantString, ":")[[1]]
                  stopifnot(length(tokens) == 4)
                  tbl.out <- data.frame(chrom=tokens[1], loc=as.numeric(tokens[2]), wt=tokens[3], mut=tokens[4])
              }
              tbl.out$chrom <- as.character(tbl.out$chrom)
              tbl.out$wt <- as.character(tbl.out$wt)
              tbl.out$mut <- as.character(tbl.out$mut)
              tbl.out
          })
#----------------------------------------------------------------------------------------------------
#' Retrieve the Sequence for a Set of Regions
#'
#' Given a MotifMatcher object, a table of chromosomal regions, and an optional set of variants,
#' return the sequences as a new column of the table.
#'
#' @rdname getSequence
#' @aliases getSequence
#'
#' @param obj An object of class MotifMatcher
#' @param tbl.regions A data frame where each row contains a chromosomal region with the fields
#' "chrom", "start", and "end".
#' @param variants A character containing variants to use for the matching (default = NA_character_)
#' The variants should either have the same number of entries as rows in the \code{tbl.regions},
#' or they should not be supplied.
#'
#' @return The \code{tbl.regions} data frame with an added column containing the sequence for each
#' entry
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Retrieve the sequences for the rs13384219 neighborhood
#' library(MotifDb)
#' motifMatcher <- MotifMatcher(genomeName="hg38",
#' pfms = as.list(query(query(MotifDb, "sapiens"), "jaspar2016")))
#' tbl.regions <- data.frame(chrom="chr2", start=57907313, end=57907333, stringsAsFactors=FALSE)
#' x <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=92)
#'
#' # Retrieve the sequences, but now include a variant
#' x.mut <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions,
#' pwmMatchMinimumAsPercentage=92, "rs13384219")
#' }

setMethod("getSequence", "MotifMatcher",

          function(obj, tbl.regions, variants=NA_character_){
              # either no variants are supplied, or a variant string is offered for each row
              stopifnot(nrow(tbl.regions) == 1)
              gr.regions <- with(tbl.regions, GRanges(seqnames=chrom, IRanges(start=start, end=end)))
              seqs <- as.character(BSgenome::getSeq(obj@genome, gr.regions))
              tbl.regions$seq <- seqs
              if(!all(is.na(variants))){  # successively inject each variant
                  tbl.variants <- do.call(rbind, lapply(variants, function(variant) .parseVariantString(obj, variant)))
                  for(rr in 1:nrow(tbl.variants))
                      tbl.regions <- .injectSnp(tbl.regions, tbl.variants[rr,])
              } # !all(is.na(variants))
              else{
                  tbl.regions$status <- rep("wt", nrow(tbl.regions))
              }
              tbl.regions
          })  # getSequence
#----------------------------------------------------------------------------------------------------
