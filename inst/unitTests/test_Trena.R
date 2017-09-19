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

   test_createGeneModel()

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

   database.filename <- system.file(package="trena", "extdata", "project.sub.db")
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
#------------------------------------------------------------------------------------------------------------------------
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

   database.filename <- system.file(package="trena", "extdata", "project.sub.db")
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
test_createGeneModel <- function()
{
   printf("--- test_createGeneModel")

   trena <- Trena("hg38")
   targetGene <- "AQP4"
   aqp4.tss <- 26865884
   chromosome <- "chr18"
   fp.source <- c("postgres://whovian/brain_hint_20")

   x <- getRegulatoryChromosomalRegions(trena, chromosome, aqp4.tss-100, aqp4.tss+100, fp.source, "AQP4", aqp4.tss)

   tbl.regulatoryRegions <- x[[fp.source]]

   tbl.geneModel <- createGeneModel(trena, targetGene, c("lasso", "randomForest"), tbl.regulatoryRegions, mtx)
   #checkTrue(all(getGeneModelTableColumnNames(trena) == colnames(tbl.geneModel)))
   tbl.strong <- subset(tbl.geneModel, rf.score > 3)
   checkTrue(all(c("TEAD1", "SP3", "KLF3", "NEUROD2") %in% tbl.strong$gene))

   checkEquals(openPostgresConnections(), 0)

   invisible(list(tbl.regulatoryRegions=tbl.regulatoryRegions,
                  tbl.geneModel=tbl.geneModel))

} # test_createGeneModel
#------------------------------------------------------------------------------------------------------------------------
# temporary hack.  the database-accessing classes should clean up after themselves
openPostgresConnections <- function()
{
   connections <- RPostgreSQL::dbListConnections(RPostgreSQL::PostgreSQL())
   length(connections)

} # openPostgresConnections
#------------------------------------------------------------------------------------------------------------------------
test_createGeneModel <- function()
{
   jaspar.human.pfms <- as.list(query(query(MotifDb, "jaspar2016"), "sapiens"))
   motifMatcher <- MotifMatcher(genomeName="hg38", pfms=jaspar.human.pfms)

      # pretend that all motifs are potentially active transcription sites - that is, ignore
      # what could be learned from open chromatin or dnasei footprints
      # use MEF2C, and 100 bases downstream, and 500 bases upstream of one of its transcripts TSS chr5:88825894

   tss <- 88825894
   tbl.region <- data.frame(chrom="chr5", start=tss-100, end=tss+400, stringsAsFactors=FALSE)
   tbl.motifs <- findMatchesByChromosomalRegion(motifMatcher, tbl.region, pwmMatchMinimumAsPercentage=92)
   tbl.motifs.tfs <- associateTranscriptionFactors(MotifDb, tbl.motifs, source="MotifDb", expand.rows=FALSE)
   solver.names <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman")
   trena <- Trena("hg38")
   tbl.geneModel <- createGeneModel(trena, "MEF2C", solver.names, tbl.motifs.tfs, mtx)

   checkTrue(is.data.frame(tbl.geneModel))

   expected.colnames <- c("gene", "beta.lasso", "lasso.p.value", "pearson.coeff", "rf.score", "beta.ridge",
                          "spearman.coeff", "concordance", "pcaMax")
   checkTrue(all(expected.colnames %in% colnames(tbl.geneModel)))
   checkTrue("FOXC1" %in% tbl.geneModel$gene)

   checkEquals(openPostgresConnections(), 0)

} # test_createGeneModel
#------------------------------------------------------------------------------------------------------------------------
