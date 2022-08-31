library(trena)
library(RUnit)
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
   test_addFeature()

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
test_addFeature <- function()
{
    message(sprintf("--- test_addFeature"))

    ft <- FeatureTable$new(target.gene="NDUFS2", reference.genome="hg38")
    ft$setFundamentalRegions(tbl.fimo)

    f <- system.file(package="trena", "extdata", "featureTable", "tbl.haploreg.european.rs4575098.RData")
    checkTrue(file.exists(f))
    tbl.haplo <- get(load(f))
    tbl.haplo$start <- tbl.haplo$hg38 - 1
    tbl.haplo$end   <- tbl.haplo$hg38

    feature.guide <- list(haploreg.rsid="rsid", haploreg.rSquared="rSquared")
    ft$addFeature(tbl.haplo, feature.genome="hg38", feature.guide)

    tbl <- ft$getTable()
    checkTrue(all(names(feature.guide) %in% colnames(tbl)))

} # test_setFundamentalRegions
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
