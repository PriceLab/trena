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
   test_addRegionFeature.haploreg()
   test_addRegionFeature.eqtls.hg19()

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

    feature.guide <- list(haploreg.rsid="rsid", haploreg.rSquared="rSquared")
    ft$addRegionFeature(tbl.haplo, feature.genome="hg38", feature.guide)

    tbl <- ft$getTable()
    checkTrue(all(names(feature.guide) %in% colnames(tbl)))

    checkTrue(nrow(subset(tbl, nchar(haploreg.rsid) > 0 & tf=="ZNF410")) >= 4)

} # test_addRegionFeature.haploreg
#----------------------------------------------------------------------------------------------------
test_addRegionFeature.eqtl.hg19 <- function()
{
    message(sprintf("--- test_addRegionFeature.eqtl.hg19"))

    tbl.eqtl <- get(load(system.file(package="trena", "extdata", "featureTable", "eqtls.ndufs2.hg19.RData")))
    tbl.eqtl <- subset(tbl.eqtl, study=="ampad-rosmap")
    dim(tbl.eqtl)  # 46 10
    checkTrue(all(c("chrom", "rsid", "hg38", "pvalue") %in% colnames(tbl.eqtl)))
    tbl.eqtl$start <- tbl.eqtl$hg38 - 1
    tbl.eqtl$end <- tbl.eqtl$hg38

    ft <- FeatureTable$new(target.gene="NDUFS2", reference.genome="hg38")
    ft$setFundamentalRegions(tbl.fimo)

    feature.guide <- list(rosmap.eqtl.rsid="rsid", rosmap.eqtl.pvalue="pvalue")
    ft$addRegionFeature(tbl.eqtl, feature.genome="hg38", feature.guide)

    tbl <- ft$getTable()
    checkTrue(all(names(feature.guide) %in% colnames(tbl)))
    table(tbl$rosmap.eqtl.rsid)
    checkTrue(nrow(subset(tbl, rosmap.eqtl.rsid=="")) > 2116200)

} # test_addRegionFeature.eqtl.hg19
#----------------------------------------------------------------------------------------------------
test_addGeneFeature.correlated.rnaseq <- function()
{
   message(sprintf("---  test_addFeature.correlated.rnaseq"))

   f <- system.file(package="trena", "extdata", "featureTable", "rna-seq.ndufs2.cor-gtex-cortex.RData")
   tbl.cor <- get(load(f))
   checkEquals(colnames(tbl.cor), c("gene", "cor"))
   checkEquals(dim(tbl.cor), c(24821, 2))

   ft <- FeatureTable$new(target.gene="NDUFS2", reference.genome="hg38")
   ft$setFundamentalRegions(tbl.fimo)

   feature.col.name <- "cortex.rna.seq.cor"
   ft$addGeneFeature(tbl.cor, feature.name=feature.col.name)

   tbl <- ft$getTable()
   checkTrue(feature.col.name %in% colnames(tbl))

      # choose a few correlations, make sure they are inscribed
   set.seed(31)
   indices <- sample(seq_len(nrow(tbl)), size=20)
   tfs.oi <- tbl$tf[indices]
   goi <- intersect(tfs.oi, tbl.cor[, "gene"])
   incoming.correlations <- tbl.cor[match(goi, tbl.cor$gene), "cor"]
   indices <- match(goi, tbl$tf)
   inscribed.correlations <- as.numeric(unlist(tbl[indices, "cortex.rna.seq.cor"]))
   checkEqualsNumeric(incoming.correlations, inscribed.correlations)
   lapply(tbl, class)

} # test_addFeature.correlated.rnaseq
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
