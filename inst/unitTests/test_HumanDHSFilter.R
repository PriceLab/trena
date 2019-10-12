library(trena)
library(MotifDb)
library(RUnit)
library(RPostgreSQL)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
mef2c.tss <- 88904257
mef2c.promoter.region <- list(chrom="chr5", start=mef2c.tss-100, end=mef2c.tss+100)
mef2c.promoter.string <- with(mef2c.promoter.region, sprintf("%s:%d-%d", chrom, start, end))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_basicConstructor()
   test_getOpenChromatinFastAndSimple()

   test_getEncodeRegulatoryTableNames()
   test_checkSampleOfEncodeTables(quiet=FALSE)

   test_getRegulatoryRegions()
   test_getCandidates.emptyRegion()
   test_getCandidates.vrk2.twoRegions()
   test_getCandidates.vrk2.rs13384219.variant()

} # runTests
#----------------------------------------------------------------------------------------------------
# reuse this in several tests
create.vrk2.candidateFilterSpec <- function(promoter.length=1000)
{
   target.gene <- "VRK2"
   genome <- "hg38"
   chromosome <- "chr2"
   tss <- 57907651
   promoter.length <- promoter.length
   db.address <- system.file(package="trena", "extdata")
   genome.db.uri    <- paste("sqlite:/", db.address, "vrk2.neighborhood.hg38.gtfAnnotation.db",  sep = "/")
   tbl.regions <- data.frame(chrom="chr2", start=57906700, end=57906870, stringsAsFactors=FALSE)

   candidateFilterSpec <- list(filterType="EncodeDNaseClusters",
                               genomeName=genome,
                               encodeTableName="wgEncodeRegDnaseClustered",
                               pwmMatchPercentageThreshold=85L,
                               geneInfoDB=genome.db.uri,
                               regions=tbl.regions,
                               pfms = as.list(query(query(MotifDb, "sapiens"),"jaspar2016")),
                               variants=NA_character_)

   candidateFilterSpec

} # create.vrk2.candidateFilterSpec
#----------------------------------------------------------------------------------------------------
# reuse this in several tests
create.vrk2.candidateFilterSpec.twoRegions <- function()
{
   target.gene <- "VRK2"
   genome <- "hg38"
   chromosome <- "chr2"
   tss <- 57907651
   promoter.length <- 1000
   db.address <- system.file(package="trena", "extdata")
   genome.db.uri    <- paste("sqlite:/", db.address, "vrk2.neighborhood.hg38.gtfAnnotation.db",  sep = "/")

   tbl.regions <- data.frame(chrom=c("chr2", "chr2"),
                             start=c(57906700, 57907740),
                             end=c(57906870, 57908150),
                             stringsAsFactors=FALSE)

   cfSpec <- list(filterType="EncodeDNaseClusters",
                  genomeName=genome,
                  encodeTableName="wgEncodeRegDnaseClustered",
                  pwmMatchPercentageThreshold=85L,
                  geneInfoDB= genome.db.uri,
                  regions=tbl.regions,
                  pfms = as.list(query(query(MotifDb, "sapiens"),"jaspar2016")),
                  variants=NA_character_)

   cfSpec

} # create.vrk2.candidateFilterSpec.twoRegions
#----------------------------------------------------------------------------------------------------
#  rs13384219  A->G
#  rs13384219 at chr2:57907073-57907573
#  gtcagtagtggtggaaccagcatgc[A/G]aattagacaatgtgacttcatagcc
#  Chromosome: 2:57907323
#  vrk2 tss at chr2:57908651
create.vrk2.rs13384219.variant.candidateFilterSpec <- function(shoulder=10)
{
   spec <- create.vrk2.candidateFilterSpec()

   rs13384219.loc <- 57907323
   left.loc <- rs13384219.loc - shoulder
   right.loc <- rs13384219.loc + shoulder
   spec$regions <- data.frame(chrom="chr2", start=left.loc, end=right.loc, stringsAsFactors=FALSE)
   spec$variants <- "rs13384219"

   spec

} # create.vrk2.rs13384219.variant.candidateFilterSpec
#----------------------------------------------------------------------------------------------------
test_basicConstructor <- function(reuse=FALSE)
{
   if(!interactive()) return(TRUE)

   if(!reuse)  printf("--- test_basicConstructor")

   candidateFilterSpec <- create.vrk2.candidateFilterSpec()

   hdf <- with(candidateFilterSpec,
               HumanDHSFilter(genomeName,
                              encodeTableName=encodeTableName,
                              pwmMatchPercentageThreshold=85L,
                              geneInfoDatabase.uri=geneInfoDB,
                              pfms = as.list(query(query(MotifDb, "sapiens"),"jaspar2016")),
                              regions=candidateFilterSpec$regions))

   checkTrue(all(c("HumanDHSFilter", "CandidateFilter") %in% is(hdf)))

      # make sure an unsupported genome triggers an error
   candidateFilterSpec$genomeName <- "intentional error in genome name"

   if(!reuse)
     checkException(hdf <- with(candidateFilterSpec,
                           HumanDHSFilter(genomeName,
                              encodeTableName=encodeTableName,
                              pwmMatchPercentageThreshold=85L,
                              geneInfoDatabase.uri=geneInfoDB,
                              pfms = as.list(query(query(MotifDb, "sapiens"),"jaspar2016")),
                              regions=regions)),
                    silent = TRUE)

   if(reuse)
      return(hdf)

} # test_basicConstructor
#----------------------------------------------------------------------------------------------------
# takes > 5 seconds, bioc check
test_getEncodeRegulatoryTableNames <- function()
{
   if(!interactive()) return(TRUE)

   printf("--- test_getEncodeRegulatoryTableNames")

   candidateFilterSpec <- create.vrk2.candidateFilterSpec()
   hdf <- with(candidateFilterSpec,
               HumanDHSFilter(genomeName,
                              encodeTableName=encodeTableName,
                              pwmMatchPercentageThreshold=85L,
                              geneInfoDatabase.uri=geneInfoDB,
                              pfms = as.list(query(query(MotifDb, "sapiens"),"jaspar2016")),
                              regions=regions))

    names <- getEncodeRegulatoryTableNames(hdf)

    checkTrue(length(names) > 90)   # 96 on (18 oct 2017)

} # test_getEncodeRegulatoryTableNames
#----------------------------------------------------------------------------------------------------
test_checkSampleOfEncodeTables <- function(quiet=TRUE)
{
   if(!interactive()) return(TRUE)

   printf("--- test_checkSampleOfEncodeTables")

   candidateFilterSpec <- create.vrk2.candidateFilterSpec()

   hdf <- with(candidateFilterSpec,
               HumanDHSFilter(genomeName,
                              encodeTableName=encodeTableName,
                              #fimoDatabase.uri=fimoDB,
                              pwmMatchPercentageThreshold=85L,
                              geneInfoDatabase.uri=geneInfoDB,
                              regions=regions,
                              pfms = as.list(query(query(MotifDb, "sapiens"),"jaspar2016")),
                              quiet=TRUE))

   tableNames <- getEncodeRegulatoryTableNames(hdf)

   #  rs13384219 at chr2:57907073-57907573
   chrom <- "chr2"
   rs13384219.loc <- 57907323
   start <- rs13384219.loc - 10000
   end <- rs13384219.loc + 10000

   selectedTableNames <- tableNames[sample(1:length(tableNames), size=10)]

   for(tableName in selectedTableNames){
      #printf("tableName: %s", tableName)
      tbl <-getRegulatoryRegions(hdf, tableName, chrom, start, end)
      if(!quiet) printf("--- %s: %d rows", tableName, nrow(tbl))
      checkTrue(nrow(tbl) >= 0)
      checkEquals(colnames(tbl), c("chrom", "chromStart", "chromEnd",  "count",  "score"))
      }

} # test_checkSampleOfEncodeTables
#----------------------------------------------------------------------------------------------------
# use this sample code to poke at the encode data offered by uscs
# note that most of the tables here only serve to list, not regions, but
# metadata:  what the inputs where, where the bb (bigBed) file can be found.
# for instance, the Hotspot table has this single line:
#                                                                     fileName
# 1 /gbdb/hg38/bbi/wgEncodeRegDnase/wgEncodeRegDnaseUwA549Hotspot.broadPeak.bb
explore.ucsc.database <- function()
{
   library(RMySQL)
   driver <- MySQL()
   host <- "genome-mysql.cse.ucsc.edu"
   user <- "genome"
   dbname <- "hg38"

   db <- dbConnect(driver, user = user, host = host, dbname = dbname)
   tables <- c("wgEncodeRegDnaseClustered", "wgEncodeRegDnaseUwA549Hotspot",  "wgEncodeRegDnaseUwA549Peak")
   main.clause <- sprintf("select * from %s where", tables[1]);

   chrom <- "chr5"
   start <- 88819630
   end   <- 88835936

   query <- paste(main.clause,
                  sprintf("chrom = '%s'", chrom),
                   sprintf("and chromStart >= %d", start),
                   sprintf("and chromEnd <= %d", end),
                   collapse = " ")
   suppressWarnings(dbGetQuery(db, sprintf("select * from %s limit 5", tables[3])))


} # explore.ucsc.database
#----------------------------------------------------------------------------------------------------
test_getRegulatoryRegions <- function()
{
   if(!interactive()) return(TRUE)

   printf("--- test_getRegulatoryRegions");

   hdf <- test_basicConstructor(reuse=TRUE)

   tableNames <- getEncodeRegulatoryTableNames(hdf)
   table <- "wgEncodeRegDnaseClustered"
   checkTrue(table %in% tableNames)


      # some DHS regions, good for testing small region overlap handling
      #    chrom chromStart chromEnd score sourceCount
      # 1   chr1   88801880 88802150   573          64
      # 2   chr1   88802480 88802930   287          13
      # 3   chr1   88803000 88803270   541          60
      # 4   chr1   88811140 88811290   100           1
      # 5   chr1   88811400 88811550    68           1

     # a small region, entirely within a DHS region
   chrom = "chr1"
   start <-  88802520
   end   <-  88802530
   tbl.regions <- getRegulatoryRegions(hdf, table, chrom, start, end)
   checkEquals(nrow(tbl.regions), 1)
   checkEquals(tbl.regions$chromStart, start)
   checkEquals(tbl.regions$chromEnd, end)

     # a small region, overhanging that DHS region at its upperbound, on the right
   chrom = "chr1"
   start <-  88802145
   end   <-  88802155
   tbl.regions <- getRegulatoryRegions(hdf, table, chrom, start, end)
   checkEquals(nrow(tbl.regions), 1)
   checkEquals(tbl.regions$chromStart, start)
   checkEquals(tbl.regions$chromEnd, 88802150)

     # a small region, overhanging that DHS region at its lowerbound, on the left
   chrom = "chr1"
   start <-  88801878
   end   <-  88801883
   tbl.regions <- getRegulatoryRegions(hdf, table, chrom, start, end)
   checkEquals(nrow(tbl.regions), 1)
   checkEquals(tbl.regions$chromStart, 88801880)
   checkEquals(tbl.regions$chromEnd, end)

      # another completely contained-in-DHS-region, small area of interest
   chrom <- "chr1"
   start <- 167830160
   end   <- 167830180
   tbl.regions <- getRegulatoryRegions(hdf, table, chrom, start, end)
   checkEquals(nrow(tbl.regions), 1)
   checkEquals(tbl.regions$chromStart, start)
   checkEquals(tbl.regions$chromEnd, end)

   # the regions around the AQP4 regulatory snps

    loc <- c(26864410, 26865469, 26855623, 26855854, 26850565)
    rsids <- c("rs3763040", "rs3875089", "rs335929", "rs3763043", "rs9951307")

    tbl.snps <- data.frame(chrom=rep("chr18", 5), start=loc-1, end=loc+1, name=rsids, stringsAsFactors=FALSE)
    tbl.regions <- getRegulatoryRegions(hdf, table, "chr18",  min(loc), max(loc))

      # with these hand-picked locations we know to expect an equal or larger number of
      # open chromatin regions in the aggregating "wgEncodeRegDnaseClustered" table
    checkTrue(nrow(tbl.regions) >= length(loc))
      # all the snp locs should fall within these regions
    checkTrue(all(loc >= tbl.regions$start))
    checkTrue(all(loc <= tbl.regions$end))

} # test_getRegulatoryRegions
#----------------------------------------------------------------------------------------------------
test_getCandidates.emptyRegion <- function()
{
   if(!interactive()) return(TRUE)

    printf("--- test_getCandidates.emptyRegion")

    db.address <- system.file(package="trena", "extdata")
    genome.db.uri <- paste("sqlite:/", db.address, "vrk2.neighborhood.hg38.gtfAnnotation.db",  sep = "/")

    hdf <-  HumanDHSFilter("hg38",
                           encodeTableName="wgEncodeRegDnaseClustered",
                           pwmMatchPercentageThreshold=80L,
                           geneInfoDatabase.uri=genome.db.uri,
                           regions=data.frame(chrom="chr18", start=26850560, end=26850565, stringsAsFactors=FALSE),
                           pfms = as.list(query(query(MotifDb, "sapiens"),"jaspar2016")),
                           quiet=TRUE)
   tbl <- getCandidates(hdf)
   checkEquals(nrow(tbl), 0)

} # test_getCandidates.emptyRegion
#----------------------------------------------------------------------------------------------------
test_getCandidates.vrk2.twoRegions <- function()
{
   if(!interactive()) return(TRUE)

   printf("--- test_getCandidates.vrk2.twoRegions")

   cfSpec <- create.vrk2.candidateFilterSpec.twoRegions()
   hdf <- with(cfSpec, HumanDHSFilter(genomeName,
                                      encodeTableName=encodeTableName,
                                      pwmMatchPercentageThreshold=97L,
                                      geneInfoDatabase.uri=geneInfoDB,
                                      region=regions,
                                      pfms = as.list(query(query(MotifDb, "sapiens"),"jaspar2016")),
                                      quiet=TRUE))

   tbl <- getCandidates(hdf)
   checkEquals(colnames(tbl),
               c("motifName", "chrom", "motifStart", "motifEnd", "strand", "motifScore", "motifRelativeScore", "match",
                 "regulatoryRegionStart", "regualtoryRegionEnd", "regulatorySequence", "variant", "shortMotif"))

    # make sure all motifs fall within the specified restions
    # cfSpec$regionsSpec: [1] "chr2:57906700-57906870" "chr2:57907740-57908150"
   starts <- tbl$motifStart
   ends   <- tbl$motifEnd
   starts.in.region.1 <- intersect(which(starts >= 57906700), which(starts <= 57906870))
   starts.in.region.2 <- intersect(which(starts >= 57907740), which(starts <= 57908150))
     # did we find them all?
   checkEquals(nrow(tbl), 11)
   checkTrue(all(sort(c(starts.in.region.1, starts.in.region.2)) == seq_len(nrow(tbl))))
   checkTrue(all(ends[starts.in.region.1] <= 57906870))
   checkTrue(all(ends[starts.in.region.2] <= 57908150))

   expected.motifs <- c("Hsapiens-jaspar2016-GATA3-MA0037.2",  "Hsapiens-jaspar2016-FOXL1-MA0033.2",
                        "Hsapiens-jaspar2016-GATA5-MA0766.1",  "Hsapiens-jaspar2016-RHOXF1-MA0719.1",
                        "Hsapiens-jaspar2016-RHOXF1-MA0719.1", "Hsapiens-jaspar2016-SPI1-MA0080.1",
                        "Hsapiens-jaspar2016-SPI1-MA0080.1",   "Hsapiens-jaspar2016-ETS1-MA0098.1",
                        "Hsapiens-jaspar2016-ETS1-MA0098.1",   "Hsapiens-jaspar2016-GATA2-MA0036.1",
                        "Hsapiens-jaspar2016-GATA2-MA0036.1")

   checkTrue(all(tbl$motifName %in% expected.motifs))

   tbl.conservative <- associateTranscriptionFactors(MotifDb, tbl, source="MotifDb", expand.rows=TRUE)
   motifDb.associated.genes <- c("GATA3","FOXL1","GATA5","RHOXF1","RHOXF1","SPI1","SPI1","ETS1","ETS1","GATA2", "GATA2")
   checkEquals(sort(motifDb.associated.genes), sort(tbl.conservative$geneSymbol))


      # we used MotifDb as our source of motifs, wherein motifName is, for example,
      # Hsapiens-jaspar2016-GATA2-MA0036.1.
      # if we use TFClass mapping from motifs to transcription factor gene names, we must create
      # a new column "shortMotif" so that the motif has a form we can use as a TFclass key
   tbl$shortMotif <- unlist(lapply(strsplit(tbl$motifName, split="-"), "[", 4))
   tbl.liberal      <- associateTranscriptionFactors(MotifDb, tbl, source="TFClass", expand.rows=TRUE)
   tfClass.associated.genes <- unique(tbl.liberal$geneSymbol)
   checkEquals(length(tfClass.associated.genes), 110)

      # despite the high number of genes associated to motifs in TFClass, not all MotifDb associations
      # are found
   checkEquals(setdiff(motifDb.associated.genes, tfClass.associated.genes), c("SPI1", "ETS1"))

} # test_getCandidates.vrk2.twoRegions
#----------------------------------------------------------------------------------------------------
test_getCandidates.vrk2.rs13384219.variant <- function()
{
   if(!interactive()) return(TRUE)

   printf("--- test_getCandidates.vrk2.rs13384219.variant")
   cfSpec <- create.vrk2.rs13384219.variant.candidateFilterSpec()

   hdf.wt <- with(cfSpec, HumanDHSFilter(genomeName,
                                         encodeTableName=encodeTableName,
                                         pwmMatchPercentageThreshold=85L,
                                         geneInfoDatabase.uri=geneInfoDB,
                                         regions=regions,
                                         variants=NA_character_,
                                         pfms = as.list(query(query(MotifDb, "sapiens"),"jaspar2016")),
                                         quiet=TRUE))

   hdf.var <- with(cfSpec, HumanDHSFilter(genomeName,
                                          encodeTableName=encodeTableName,
                                          pwmMatchPercentageThreshold=85L,
                                          geneInfoDatabase.uri=geneInfoDB,
                                          regions=regions,
                                          variants=variants,
                                          pfms = as.list(query(query(MotifDb, "sapiens"),"jaspar2016")),
                                          quiet=TRUE))

   tbl.wt <- getCandidates(hdf.wt)
   tbl.var <- getCandidates(hdf.var)

   checkEquals(dim(tbl.wt),  c(53, 13))
   checkEquals(dim(tbl.var), c(45, 13))

      # many motifs survive, 9 are lost, 2 are gained
      # see trena::assessSnp for close comparison of losses and gains, wherein differing scores are reported

   checkEquals(length(intersect(tbl.wt$motifName, tbl.var$motifName)), 36)
   checkEquals(length(setdiff(tbl.wt$motifName, tbl.var$motifName)),    9)
   checkEquals(length(setdiff(tbl.var$motifName, tbl.wt$motifName)),    2)

} # test_getCandidates.vrk2.rs13384219.variant
#------------------------------------------------------------------------------------------------------------------------
# demonstrate an initially unanticipated use of this class: just get the open chromatin regions from ucsc.  a hack...
test_getOpenChromatinFastAndSimple <- function()
{
   printf("--- test_getOpenChromatinFastAndSimple")

   hdf <- HumanDHSFilter("hg38", "wgEncodeRegDnaseClustered", pwmMatchPercentageThreshold=0,
                         geneInfoDatabase.uri="bogus", regions=data.frame(), pfms=list())
   chrom <- "chr14"
   start <- 91654248
   end <-   93010778
   tbl.dhs <- getRegulatoryRegions(hdf, "wgEncodeRegDnaseClustered", chrom, start, end)

   checkEquals(colnames(tbl.dhs), c("chrom", "chromStart", "chromEnd", "count", "score"))
   checkTrue(nrow(tbl.dhs) > 1200)
   checkTrue(all(tbl.dhs$chrom == chrom))
   checkTrue(all(tbl.dhs$start >= start))
   checkTrue(all(tbl.dhs$end <= end))

} # test_getOpenChromatinFastAndSimple
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
