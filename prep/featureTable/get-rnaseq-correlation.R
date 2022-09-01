library(EndophenotypeExplorer)
targetGene <- "NDUFS2"
etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD")
tbl.eqtl <- etx$get.ampad.EQTLsForGene()
tbl.eqtl.ndufs2.hg19 <- subset(tbl.eqtl, pvalue < 1e-4)
save(tbl.eqtl.ndufs2.hg19, file="../../inst/extdata/featureTable/eqtls.ndufs2.hg19.RData")
