library(trena)
library(RUnit)
library(ghdb)
#----------------------------------------------------------------------------------------------------
f <- system.file(package="trena", "extdata", "featureTable", "tbl.fimo-ndufs2-47kb.RData")
tbl.fimo <- get(load(f))
tbl.fimo <- get(load("~/github/tms-makeAndBreak/studies/rs4575098-ndufs2-region/shared/tbl.fimo.NDUFS2.RData"))
p.value.col <- grep("^p\\.value$", colnames(tbl.fimo))
if(length(p.value.col) == 1)
    colnames(tbl.fimo)[p.value.col] <- "fimo.pval"
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
   test_setFundamentalRegions()
   test_addRegionFeature.haploreg()
   test_addRegionFeature.eqtls.hg19()
   test_addGeneFeature.correlated.rnaseq()
   test_addGeneHancer()
   test_addDistanceToTSS()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    ft <- FeatureTable$new(target.gene="NDUFS2", reference.genome="hg38")
    checkTrue(all(c("FeatureTable", "R6") %in% class(ft)))

    tbl <- ft$getTable()
    checkEquals(dim(tbl), c(0,0))

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_setFundamentalRegions <- function()
{
    message(sprintf("--- test_eetFundamentalRegions"))

    ft <- FeatureTable$new(target.gene="NDUFS2", reference.genome="hg38")

    ft$setFundamentalRegions(tbl.fimo)
    #ft$setFundamentalRegions(tbl.fimo)
    tbl <- ft$getTable()
    checkTrue(ncol(tbl) == 9)
    checkTrue(nrow(tbl) > 2000000)


} # test_setFundamentalRegions
#----------------------------------------------------------------------------------------------------
test_addRegionFeature.haploreg <- function()
{
    message(sprintf("--- test_addRegionFeature.haploreg"))

    ft <- FeatureTable$new(target.gene="NDUFS2", reference.genome="hg38")
    ft$setFundamentalRegions(tbl.fimo)

    f <- system.file(package="trena", "extdata", "featureTable", "tbl.haploreg.european.rs4575098.RData")
    checkTrue(file.exists(f))
    tbl.haplo <- get(load(f))
    tbl.haplo$start <- tbl.haplo$hg38 - 1
    tbl.haplo$end   <- tbl.haplo$hg38

    feature.guide <-  list(haploreg.rsid="rsid", haploreg.rSquared="rSquared")
    default.values <- list(haploreg.rsid="",     haploreg.rSquared=0)
    ft$addRegionFeature(tbl.haplo, feature.genome="hg38", feature.guide, default.values)

    tbl <- ft$getTable()
    checkTrue(all(names(feature.guide) %in% colnames(tbl)))

    checkTrue(nrow(subset(tbl, nchar(haploreg.rsid) > 0 & tf=="ZNF410")) >= 4)

} # test_addRegionFeature.haploreg
#----------------------------------------------------------------------------------------------------
test_addRegionFeature.eqtls.hg19 <- function()
{
    message(sprintf("--- test_addRegionFeature.eqtls.hg19"))

    tbl.eqtl <- get(load(system.file(package="trena", "extdata", "featureTable", "eqtls.ndufs2.hg19.RData")))
    tbl.eqtl <- subset(tbl.eqtl, study=="ampad-rosmap")
    dim(tbl.eqtl)  # 46 10
    checkTrue(all(c("chrom", "rsid", "hg38", "pvalue") %in% colnames(tbl.eqtl)))
    tbl.eqtl$start <- tbl.eqtl$hg38 - 1
    tbl.eqtl$end <- tbl.eqtl$hg38

    ft <- FeatureTable$new(target.gene="NDUFS2", reference.genome="hg38")
    ft$setFundamentalRegions(tbl.fimo)

    feature.guide <- list(rosmap.eqtl.rsid="rsid", rosmap.eqtl.pvalue="pvalue")
    default.values <- list(rosmap.eqtl.rsid="", rosmap.eqtl.pvalue=1)
    ft$addRegionFeature(tbl.eqtl, feature.genome="hg38", feature.guide, default.values)

    tbl <- ft$getTable()
    checkTrue(all(names(feature.guide) %in% colnames(tbl)))
    table(tbl$rosmap.eqtl.rsid)
    checkTrue(nrow(subset(tbl, rosmap.eqtl.rsid=="")) > 2116200)

} # test_addRegionFeature.eqtls.hg19
#----------------------------------------------------------------------------------------------------
test_addGeneFeature.correlated.rnaseq <- function()
{
   message(sprintf("--- test_addFeature.correlated.rnaseq"))

   f <- system.file(package="trena", "extdata", "featureTable", "rna-seq.ndufs2.cor-gtex-cortex.RData")
   tbl.cor <- get(load(f))
   checkEquals(colnames(tbl.cor), c("gene", "cor"))
   checkEquals(dim(tbl.cor), c(24821, 2))

   ft <- FeatureTable$new(target.gene="NDUFS2", reference.genome="hg38")
   ft$setFundamentalRegions(tbl.fimo)

   feature.col.name <- "gtex.cortex.rna.seq.cor"
   ft$addGeneFeature(tbl.cor, feature.name=feature.col.name, default.value=0)

   tbl <- ft$getTable()
   checkTrue(feature.col.name %in% colnames(tbl))

      # choose a few correlations, make sure they are inscribed
   set.seed(31)
   indices <- sample(seq_len(nrow(tbl)), size=20)
   tfs.oi <- tbl$tf[indices]
   goi <- intersect(tfs.oi, tbl.cor[, "gene"])
   incoming.correlations <- tbl.cor[match(goi, tbl.cor$gene), "cor"]
   indices <- match(goi, tbl$tf)
   inscribed.correlations <- as.numeric(unlist(tbl[indices, "gtex.cortex.rna.seq.cor"]))
   checkEqualsNumeric(incoming.correlations, inscribed.correlations)
   lapply(tbl, class)

} # test_addFeature.correlated.rnaseq
#----------------------------------------------------------------------------------------------------
test_addGeneHancer <- function()
{
    message(sprintf("--- test_addGeneHancer"))

    ghdb <- GeneHancerDB()
    goi <- "NDUFS2"
    tbl.gh <- retrieveEnhancersFromDatabase(ghdb, goi, tissues="brain")

    ft <- FeatureTable$new(target.gene="NDUFS2", reference.genome="hg38")
    ft$setFundamentalRegions(tbl.fimo)

    feature.guide <- list(gh.brain="tissue", gh.brain.score="combinedscore")
    default.values <- list(gh.brain="", gh.brain.score=0)
    feature.guide <- list(gh.brain.score="combinedscore")
    default.values <- list(gh.brain.score=0)

    ft$addRegionFeature(tbl.gh, feature.genome="hg38", feature.guide, default.values)

    tbl <- ft$getTable()
    checkEqualsNumeric(fivenum(tbl$gh.brain.score), c(0, 0, 0, 0, 261), tol=1)
    checkEquals(nrow(subset(tbl, gh.brain.score > 260)), 28220)
    checkEquals(nrow(subset(tbl, gh.brain.score > 10 & gh.brain.score < 20)), 0)

    tbl.gh <- retrieveEnhancersFromDatabase(ghdb, goi, tissues="all")
    tbl.gh$tissue <- substring(tbl.gh$tissue, 1, 10)
    feature.guide <-  list(gh.any.tissue="tissue", gh.any.tissue.score="combinedscore")
    default.values <- list(gh.any.tissue="",       gh.any.tissue.score=0)
    ft$addRegionFeature(tbl.gh, feature.genome="hg38", feature.guide, default.values)

    tbl <- ft$getTable()
    checkEquals(nrow(subset(tbl, gh.any.tissue.score > 260)), 28220)
    checkTrue(nrow(subset(tbl, gh.any.tissue.score > 10 & gh.any.tissue.score < 20)) > 30000)

} # test_addGeneHancer
#----------------------------------------------------------------------------------------------------
test_addDistanceToTSS <- function()
{
   message(sprintf("--- test_addDistanceToTSS"))

   ft <- FeatureTable$new(target.gene="NDUFS2", reference.genome="hg38")
   ft$setFundamentalRegions(tbl.fimo)

   ft$addDistanceToTSS(161197104)
   tbl <- ft$getTable()
   checkEqualsNumeric(fivenum(tbl$distanceToTSS), c(-512978, -243587, -7591, 243276, 491713), tol=100)
   checkTrue(nrow(subset(tbl, distanceToTSS == 0)) > 5)   # 8 seen on (1 sep 2022)

} # test_addDistanceToTSS
#----------------------------------------------------------------------------------------------------
demo <- function()
{
   ft <- FeatureTable$new(target.gene="NDUFS2", reference.genome="hg38")
   ft$setFundamentalRegions(tbl.fimo)

   f <- system.file(package="trena", "extdata", "featureTable", "tbl.haploreg.european.rs4575098.RData")
   checkTrue(file.exists(f))
   tbl.haplo <- get(load(f))
   tbl.haplo$start <- tbl.haplo$hg38 - 1
   tbl.haplo$end   <- tbl.haplo$hg38

   feature.guide <- list(haploreg.rsid="rsid", haploreg.rSquared="rSquared")
   ft$addRegionFeature(tbl.haplo, feature.genome="hg38", feature.guide)

   tbl.eqtl <- get(load(system.file(package="trena", "extdata", "featureTable", "eqtls.ndufs2.hg19.RData")))
   tbl.eqtl <- subset(tbl.eqtl, study=="ampad-rosmap")
   dim(tbl.eqtl)  # 46 10
   checkTrue(all(c("chrom", "rsid", "hg38", "pvalue") %in% colnames(tbl.eqtl)))
   tbl.eqtl$start <- tbl.eqtl$hg38 - 1
   tbl.eqtl$end <- tbl.eqtl$hg38

   feature.guide <- list(rosmap.eqtl.rsid="rsid", rosmap.eqtl.pvalue="pvalue")
   ft$addRegionFeature(tbl.eqtl, feature.genome="hg38", feature.guide)

   f <- system.file(package="trena", "extdata", "featureTable", "rna-seq.ndufs2.cor-gtex-cortex.RData")
   tbl.cor <- get(load(f))
   checkEquals(colnames(tbl.cor), c("gene", "cor"))
   checkEquals(dim(tbl.cor), c(24821, 2))

   feature.col.name <- "cortex.rna.seq.cor"
   ft$addGeneFeature(tbl.cor, feature.name=feature.col.name)

   tbl <- ft$getTable()
   dim(tbl)
   colnames(tbl)
   coi <- c("tf", "cortex.rna.seq.dor", "rosmap.eqtl.pvalue", "chrom", "start", "end")
   tbl <- tbl[, ..coi]
   tbl$rosmap.eqtl.pvalue <- -log10(tbl$rosmap.eqtl.pvalue+000001)
   colnames(tbl)[6] <- "rosmap.eqtl.score"
   dim(tbl)

   save(tbl, file="tbl.RData")

} # demo
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
