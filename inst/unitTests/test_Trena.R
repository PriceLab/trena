library(trena)
library(RUnit)
library(MotifDb)
library(RPostgreSQL)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
if(!exists("mtx")){
    load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
    mtx <- asinh(mtx.sub)
    mtx.var <- apply(mtx, 1, var)
    deleters <- which(mtx.var < 0.01)
    if(length(deleters) > 0)   # 15838 x 638
        mtx <- mtx[-deleters,]
}
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_basicConstructor()

    test_getRegulatoryRegions_oneFootprintSource()
    test_getRegulatoryRegions_encodeDHS()
    test_getRegulatoryRegions_twoFootprintSources()

    test_createGeneModelFromRegulatoryRegions()
    test_createGeneModelFromTfList()

    test_getProximalPromoterHuman()
    test_getProximalPromoterMouse()

    checkEquals(openPostgresConnections(), 0)

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_basicConstructor <- function()
{
    printf("--- test_basicConstructor")

    trena <- Trena("hg38")
    checkEquals(is(trena), "Trena")

    checkException(Trena("hg00"), silent=TRUE)
    checkException(Trena(""), silent=TRUE)
    checkEquals(openPostgresConnections(), 0)

} # test_basicConstructor
#----------------------------------------------------------------------------------------------------
test_getRegulatoryRegions_oneFootprintSource <- function()
{
    printf("--- test_getRegulatoryRegions_oneFootprintSource")

    trena <- Trena("hg38")
    # the package's demo sqlite database is limited to in and around hg38 MEF2C
    chromosome <- "chr5"
    mef2c.tss <- 88904257   # minus strand

    database.filename <- system.file(package="trena", "extdata", "mef2c.neigborhood.hg38.footprints.db")
    database.uri <- sprintf("sqlite://%s", database.filename)
    sources <- c(database.uri)

    loc.start <- mef2c.tss - 1000
    loc.end   <- mef2c.tss + 1000

    x <- getRegulatoryChromosomalRegions(trena, chromosome, mef2c.tss-1000, mef2c.tss+1000, sources, "MEF2C", mef2c.tss)

    checkTrue(is(x, "list"))
    # only one source, thus only one element in the result list
    checkEquals(length(x), 1)
    checkEquals(names(x), as.character(sources))
    checkTrue(all(unlist(lapply(x, function(element) is(element, "data.frame")), use.names=FALSE)))

    tbl.reg <- x[[sources[[1]]]]
    checkEquals(colnames(tbl.reg), getRegulatoryTableColumnNames(trena))

    checkEquals(unique(tbl.reg$chrom), chromosome)
    checkTrue(nrow(tbl.reg) > 20)

    checkTrue(all(tbl.reg$motifStart >= loc.start))
    checkTrue(all(tbl.reg$motifStart <= loc.end))
    checkTrue(all(tbl.reg$motifEnd   >= loc.start))
    checkTrue(all(tbl.reg$motifEnd   <= loc.end))
    checkTrue(all(tbl.reg$distance.from.tss >= -1000))
    checkTrue(all(tbl.reg$distance.from.tss <=  1000))

    checkEquals(openPostgresConnections(), 0)

} # test_getRegulatoryRegions_oneFootprintSource
#----------------------------------------------------------------------------------------------------
test_getRegulatoryRegions_encodeDHS <- function()
{
    printf("--- test_getRegulatoryRegions_encodeDHS")

    trena <- Trena("hg38")
    #chromosome <- "chr5"
    #mef2c.tss <- 88904257   # minus strand
    #loc.start <- mef2c.tss - 10000
    #loc.end   <- mef2c.tss + 10000
    #sources <- c("encodeHumanDHS")

    aqp4.tss <- 26865884
    chromosome <- "chr18"
    sources <- c("encodeHumanDHS")

    # first submit a tiny region with no regulatory regions
    tbl <- getRegulatoryChromosomalRegions(trena, chromosome, aqp4.tss-1, aqp4.tss+3, sources, "AQP4", aqp4.tss)
    checkEquals(nrow(tbl[["encodeHumanDHS"]]), 0)

    # now a larger region
    x <- getRegulatoryChromosomalRegions(trena, chromosome, aqp4.tss-100, aqp4.tss+100, sources, "AQP4", aqp4.tss)

    checkTrue(is(x, "list"))
    checkEquals(names(x), as.character(sources))
    checkTrue(all(unlist(lapply(x, function(element) is(element, "data.frame")), use.names=FALSE)))

    tbl.reg <- x[[sources[[1]]]]
    checkTrue(all(colnames(tbl.reg) == getRegulatoryTableColumnNames(trena)))
    checkTrue(nrow(tbl.reg) > 20)
    checkEquals(length(grep("AQP4.dhs", tbl.reg$id)), nrow(tbl.reg))

    checkEquals(openPostgresConnections(), 0)

} # test_getRegulatoryRegions_encodeDHS
#----------------------------------------------------------------------------------------------------
# for quick testing, we use only the small (mef2c-centered) sqlite database distributed
# with the package.  so this "two-source" test uses that one source twice, producing
# a list of length 2, each with the same name, each with the same contents
test_getRegulatoryRegions_twoFootprintSources <- function()
{
    printf("--- test_getRegulatoryRegions_twoFootprintSources")

    trena <- Trena("hg38")
    # the package's demo sqlite database is limited to in and around hg38 MEF2C
    chromosome <- "chr5"
    mef2c.tss <- 88904257   # minus strand

    database.filename <- system.file(package="trena", "extdata", "mef2c.neigborhood.hg38.footprints.db")
    database.uri <- sprintf("sqlite://%s", database.filename)
    sources <- c(database.uri,  database.uri)

    loc.start <- mef2c.tss - 1000
    loc.end   <- mef2c.tss + 1000

    x <- getRegulatoryChromosomalRegions(trena, chromosome, mef2c.tss-1000, mef2c.tss+1000, sources, "MEF2C", mef2c.tss)

    checkTrue(is(x, "list"))
    checkEquals(length(x), 2)
    checkEquals(names(x), sources)
    checkTrue(all(unlist(lapply(x, function(element) is(element, "data.frame")), use.names=FALSE)))

    tbl.reg <- x[[sources[[1]]]]
    checkEquals(colnames(tbl.reg), getRegulatoryTableColumnNames(trena))

    checkEquals(unique(tbl.reg$chrom), chromosome)
    checkTrue(nrow(tbl.reg) > 20)

    checkTrue(all(tbl.reg$motifStart >= loc.start))
    checkTrue(all(tbl.reg$motifStart <= loc.end))
    checkTrue(all(tbl.reg$motifEnd   >= loc.start))
    checkTrue(all(tbl.reg$motifEnd   <= loc.end))
    checkTrue(all(tbl.reg$distance.from.tss >= -1000))
    checkTrue(all(tbl.reg$distance.from.tss <=  1000))

    checkEquals(openPostgresConnections(), 0)

} # test_getRegulatoryRegions_twoFootprintSources
#----------------------------------------------------------------------------------------------------
# temporary hack.  the database-accessing classes should clean up after themselves
openPostgresConnections <- function()
{
    connections <- RPostgreSQL::dbListConnections(RPostgreSQL::PostgreSQL())
    length(connections)

} # openPostgresConnections
#----------------------------------------------------------------------------------------------------
test_createGeneModelFromRegulatoryRegions <- function()
{
    printf("--- test_createGeneModelFromRegulatoryRegions")

    targetGene <- "MEF2C"
    jaspar.human.pfms <- as.list(query(MotifDb, "jaspar2016", "sapiens"))
    motifMatcher <- MotifMatcher(genomeName="hg38", pfms=jaspar.human.pfms)

    # pretend that all motifs are potentially active transcription sites - that is, ignore
    # what could be learned from open chromatin or dnasei footprints
    # use MEF2C, and 100 bases downstream, and 500 bases upstream of one of its transcripts TSS chr5:88825894

    tss <- 88825894
    tbl.region <- data.frame(chrom="chr5", start=tss-100, end=tss+500, stringsAsFactors=FALSE)
    tbl.region <- data.frame(chrom="chr5", start=tss-5000, end=tss+5000, stringsAsFactors=FALSE)
    tbl.motifs <- findMatchesByChromosomalRegion(motifMatcher, tbl.region, pwmMatchMinimumAsPercentage=85)
    tbl.motifs.tfs <- associateTranscriptionFactors(MotifDb, tbl.motifs, source=c("MotifDb", "TFClass"), expand.rows=FALSE)
    solver.names <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman")
    trena <- Trena("hg38")

       # very meager model results - but useful for testing nonetheless.
    tbl.geneModel <- createGeneModelFromRegulatoryRegions(trena, targetGene, solver.names, tbl.motifs.tfs, mtx)

    checkTrue(is.data.frame(tbl.geneModel))

    expected.colnames <- c("gene", "betaLasso", "lassoPValue", "pearsonCoeff", "rfScore", "betaRidge",
                           "spearmanCoeff", "bindingSites")
    checkTrue(all(expected.colnames %in% colnames(tbl.geneModel)))
    checkTrue(nrow(tbl.geneModel) > 40)
    checkTrue("HLF" %in% tbl.geneModel$gene[1:10])

    checkEquals(openPostgresConnections(), 0)

} # test_createGeneModelFromRegulatoryRegions
#----------------------------------------------------------------------------------------------------
test_createGeneModelFromTfList <- function()
{
    printf("--- test_createGeneModelFromTfList")

    targetGene <- "MEF2C"

    solver.names <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman")

    candidate.tfs <-   c("HLF", "STAT4", "SATB2", "SATB1", "TSHZ3", "TSHZ2", "FOXP2",
                         "FOXP1", "LHX6", "BACH1", "SOX12", "FOXD4L1", "NFE2L2", "ZHX3",
                         "ZBTB16", "ZHX1", "TAF1", "STAT6", "POU4F1", "FOXD1", "ATF2",
                         "BCL6B", "STAT5B", "NR5A2", "FOXE3", "STAT3", "ATF7", "STAT2")

    trena <- Trena("hg38")

    tbl.geneModel <- createGeneModelFromTfList(trena, targetGene, solver.names, candidate.tfs, mtx)
    checkTrue(is.data.frame(tbl.geneModel))

    expected.colnames <- c("gene", "betaLasso", "lassoPValue", "pearsonCoeff", "rfScore", "betaRidge",
                           "spearmanCoeff",  "bindingSites")
    checkTrue(all(expected.colnames %in% colnames(tbl.geneModel)))
    checkTrue(nrow(tbl.geneModel) > 25)
    checkTrue(all(c("HLF", "STAT4", "SATB2") %in% tbl.geneModel$gene[1:10]))

    checkEquals(openPostgresConnections(), 0)

} # test_createGeneModelFromTfList
#----------------------------------------------------------------------------------------------------
test_getProximalPromoterHuman <- function()
{
    printf("--- test_getProximalPromoterHuman")

    if(!interactive())
       return(TRUE)

    trena <- Trena("hg38")

    # Designate the MEF2C gene and a shoulder size of 1000
    geneSymbol <- "MEF2C"
    tssUpstream <- 1000
    tssDownstream <- 1000

       # Pull the regions for MEF2C
    regions <- getProximalPromoter(trena, geneSymbol, tssUpstream, tssDownstream)

       # Check the type of data returned and its size
    checkEquals(dim(regions), c(1,4))
    checkEquals(class(regions), "data.frame")
    checkEquals(with(regions,  end - start), tssUpstream + tssDownstream)

      # check with bogus gene symbol
    checkTrue(is.na(getProximalPromoter(trena, "bogus", tssUpstream, tssDownstream)))

    # Test it on a list, with one real and one bogus, and make sure its the same as the first
    genes.of.interest <- c(geneSymbol, "bogus")
    regions.2 <- getProximalPromoter(trena, genes.of.interest, tssUpstream, tssDownstream)
    checkEquals(regions, regions.2)

} # test_getProximalPromoterHuman
#----------------------------------------------------------------------------------------------------
test_getProximalPromoterMouse <- function(){

    printf("--- test_getProximalPromoterMouse")

    if(!interactive())
       return(TRUE)

    trena <- Trena("mm10")

    # Designate the Twist2 gene and a shoulder size of 1000
    geneSymbol <- "Twist2"
    tssUpstream <- 1000
    tssDownstream <- 1000

    # Pull the regions for Twist2
    regions <- getProximalPromoter(trena, geneSymbol, tssUpstream, tssDownstream)

    # Check the type of data returned and its size
    checkEquals(dim(regions), c(1,4))
    checkEquals(class(regions), "data.frame")

    # Check the nominal values (tss = 88904257)
    tss <- 91801461
    checkEquals(regions$chrom, "chr1")
    checkEquals(regions$start, tss - tssUpstream)
    checkEquals(regions$end, tss + tssDownstream)

    # check with bogus gene symbol
    checkTrue(is.na(getProximalPromoter(trena, "bogus", tssUpstream, tssDownstream)))

    # Test it on a list, with one real and one bogus, and make sure its the same as the first
    regions.2 <- getProximalPromoter(trena, c(geneSymbol, "bogus"), tssUpstream, tssDownstream)
    checkEquals(regions, regions.2)

} # test_getProximalPromoterMouse
#----------------------------------------------------------------------------------------------------
test_assessSnp <- function()
{
   if(!interactive()) return()

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
#----------------------------------------------------------------------------------------------------
# in preparation for adding, and ongoing testing of, a delta column for all entries, here we use a snp which at the 80%
# match level # returns all three kinds of match: in.both. wt.only, mut.only
test_assessSnp_allTypesWithDeltas <- function()
{
    if(!interactive()) return()

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
if(!interactive()) runTests()
