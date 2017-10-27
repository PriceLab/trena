setGeneric('assessSnp', signature='obj', function(obj, pfms, variant, shoulder, pwmMatchMinimumAsPercentage, relaxedMatchDelta=25)
    standardGeneric('assessSnp'))
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
#'

setMethod('assessSnp', 'Trena',

          function(obj, pfms, variant, shoulder, pwmMatchMinimumAsPercentage, relaxedMatchDelta=25){

              motifMatcher <- MotifMatcher(genomeName=obj@genomeName, pfms=pfms, quiet=obj@quiet)
              tbl.variant <- try(.parseVariantString(motifMatcher, variant), silent=TRUE)
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
                                                        variant=variant)
              if(nrow(tbl.mut) == 0){
                  warning(sprintf("no motifs altered by %s with shoulder %d", variant, shoulder))
                  return(data.frame())
              }

              tbl.wt$signature <- sprintf("%s;%s;%s", tbl.wt$motifName, tbl.wt$motifStart, tbl.wt$strand)
              tbl.mut$signature <- sprintf("%s;%s;%s", tbl.mut$motifName, tbl.mut$motifStart, tbl.mut$strand)

              # comine wt and mut tables, reorder columns and rows for easier comprehension
              tbl <- rbind(tbl.wt[, c(1,12,2,3,4,5,6,7,8,13)], tbl.mut[, c(1,12,2,3,4,5,6,7,8,13)])
              tbl <- tbl[order(tbl$motifName, tbl$motifRelativeScore, decreasing=TRUE),]
              #tbl$signature <- sprintf("%s;%s;%s", tbl$motifName, tbl$motifStart, tbl$strand)
              tbl <- tbl[,c(1,2,3:10)]

              # now look for less stringent matches.  these will be matched up with the
              # wt and mut motifs which do not yet have partners, thus enabling us to
              # provide a wt->mut motifScore.delta for each
              relaxedMatchPercentage <- pwmMatchMinimumAsPercentage-relaxedMatchDelta
              tbl.wt.relaxed <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, relaxedMatchPercentage)
              tbl.wt.relaxed$signature <- sprintf("%s;%s;%s", tbl.wt.relaxed$motifName, tbl.wt.relaxed$motifStart, tbl.wt.relaxed$strand)
              tbl.mut.relaxed <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, relaxedMatchPercentage, variant=variant)
              tbl.mut.relaxed$signature <- sprintf("%s;%s;%s", tbl.mut.relaxed$motifName, tbl.mut.relaxed$motifStart, tbl.mut.relaxed$strand)

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
# Unit tests and such

test_assessSnp <- function()
{
    printf("--- test_assessSnp")

    trena <- Trena("hg38")
    jaspar.human.pfms <- as.list(query(query(MotifDb, "jaspar2016"), "sapiens"))

    # first check for bogus variant name
    bogus.variant <- "rsBogus"
    checkEquals(assessSnp(trena, jaspar.human.pfms, bogus.variant, shoulder=5, pwmMatchMinimumAsPercentage=65),
                data.frame())

    variant <- "rs3875089"   # chr18:26865469  T->C

    # a shoulder of 3 gives us a search region of chr18:26865466-26865472
    shoulder <- 3
    # a 65% match is relaxed enough to get these results, good fodder for testing. tbl.wt then tbl.mut

    #                           motifName chrom motifStart motifEnd strand motifScore motifRelativeScore
    #   Hsapiens-jaspar2016-ETS1-MA0098.1 chr18   26865467 26865472      +   3.825000          0.8843931
    #   Hsapiens-jaspar2016-SPI1-MA0080.1 chr18   26865467 26865472      -   3.684211          0.7865169
    #  Hsapiens-jaspar2016-GATA2-MA0036.1 chr18   26865467 26865471      -   3.509434          0.9489796
    #  Hsapiens-jaspar2016-GATA3-MA0037.1 chr18   26865466 26865471      -   3.047619          0.6808511
    #  Hsapiens-jaspar2016-GATA2-MA0036.1 chr18   26865466 26865470      +   2.547170          0.6887755
    #
    #                            motifName chrom motifStart motifEnd strand motifScore motifRelativeScore
    # Hsapiens-jaspar2016-ZNF354C-MA0130.1 chr18   26865467 26865472      +    3.50000          0.6913580
    #    Hsapiens-jaspar2016-ETS1-MA0098.1 chr18   26865467 26865472      +    2.87500          0.6647399
    #   Hsapiens-jaspar2016-GATA2-MA0036.1 chr18   26865467 26865471      -    2.54717          0.6887755
    #
    # ma0098+  wt=0.884  mut=0.664
    # ma0036-  wt=0.948 only
    # ma0036+  wt=mut=0.688
    # ma0037-  wt only
    # ma0130+  mut only
    # ma0080-  wt only

    tbl.assay <- assessSnp(trena, jaspar.human.pfms, variant, shoulder, pwmMatchMinimumAsPercentage=65)
    checkEquals(dim(tbl.assay), c(8, 12))

    expected.colnames <- c("motifName", "status", "assessed", "motifRelativeScore", "delta", "signature",
                           "chrom", "motifStart", "motifEnd", "strand", "match", "variant")

    checkEquals(sort(colnames(tbl.assay)), sort(expected.colnames))
    # pull out crucial columns for checking
    tbl.test <- tbl.assay[, c("signature", "status", "assessed", "motifRelativeScore", "delta")]

    # all 3 categories should be present
    checkEquals(as.list(table(tbl.test$assessed)), list(in.both=4, mut.only=1, wt.only=3))

    # deltas are zero if both wt and mut for a motif/strand were found
    # in this case, the delta can be read off the two motifRelativeScore valus
    checkTrue(all(tbl.test$delta[grep("both", tbl.test$assessed)]==0))
    checkTrue(all(tbl.test$delta[grep("only", tbl.test$assessed)]!=0))

    tbl.ma0098 <- tbl.test[grep("MA0098.1;26865467;+", tbl.test$signature, fixed=TRUE),]
    checkEquals(nrow(tbl.ma0098), 2)
    checkTrue(all(tbl.ma0098$delta == 0))
    checkTrue(all(c("wt", "mut") %in% tbl.ma0098$status))
    checkEqualsNumeric(tbl.test[grep("MA0098", tbl.test$signature),"motifRelativeScore"],
                       c(0.8843931, 0.6647399), tol=1e-5)

    # now test for an empty table - no wt or mut motifs for this region at this minimum match
    suppressWarnings(tbl.assay.short <- assessSnp(trena, jaspar.human.pfms, "rs3875089", 3, pwmMatchMinimumAsPercentage=95))
    checkEquals(nrow(tbl.assay.short), 0)

} # test_assessSnp
#------------------------------------------------------------------------------------------------------------------------
# in preparation for adding, and ongoing testing of, a delta column for all entries, here we use a snp which at the 80%
# match level # returns all three kinds of match: in.both. wt.only, mut.only
test_assessSnp_allTypesWithDeltas <- function()
{
    printf("--- test_assessSnp_allTypesWithDeltas")

    trena <- Trena("hg38")
    jaspar.human.pfms <- as.list(query(query(MotifDb, "jaspar2016"), "sapiens"))
    snp <- "rs3763043"
    shoulder <- 10

    tbl.assay <- assessSnp(trena, jaspar.human.pfms, snp, shoulder=shoulder, pwmMatchMinimumAsPercentage=80)
    checkEquals(sort(unique(tbl.assay$assessed)), c("in.both", "mut.only", "wt.only"))

    checkEquals(ncol(tbl.assay), 12)
    checkTrue("delta" %in% colnames(tbl.assay))
    checkEqualsNumeric(min(tbl.assay$delta), -0.0487, tol=1e-3)
    checkEqualsNumeric(max(tbl.assay$delta),  0.22, tol=1e-2)

} # test_assessSnp_allTypesWithDeltas
#----------------------------------------------------------------------------------------------------
# Extra bit from MotifMatcher

