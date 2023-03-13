library(trena)
library(RUnit)
targetGene <- "ZBTB7A"
#----------------------------------------------------------------------------------------------------
fetch.eqtls <- function()
{
   require(ebi.eqtls)
   ee <- ebi.eqtls$new()
   tbl.cat <- ee$getCatalog()
   toi <- subset(tbl.cat, Tissue_group %in% c("blood", "Whole Blood") &
                          quant_method=="ge")$unique_id
      #  "Lepik_2017.blood" "TwinsUK.blood"    "GTEx.blood" "GTEx_V8.Whole_Blood"

   start.144k <- 3987816
   end.144k <-   4132224
   start.long <- 2652636
   end.long   <- 4719616

   chrom <- "chr19"
   loc.start <- start.144k
   loc.end   <- end.144k

   loc.start <- start.long
   loc.end <- end.long

   sprintf("%5.2fk", round((loc.end - loc.start)/1000, digits=2))
   fetch <- function(tissue){
       ee$fetch.eqtls.in.chunks(chrom="chr19",
                                start=loc.start,
                                end=loc.end,
                                study=tissue,
                                simplify=TRUE,
                                chunk.size=5000)
       } # fetch
    tbls <- lapply(toi, fetch)
    tbl.eqtl.gtex <- do.call(rbind, tbls)
    rownames(tbl.eqtl.gtex) <- NULL
    dim(tbl.eqtl.gtex) # 94472
    head(tbl.eqtl.gtex)

    tbl.eqtl.gtex <- subset(tbl.eqtl.gtex, gene==targetGene & pvalue < 0.2)
    dim(tbl.eqtl.gtex)  # 602
    fivenum(tbl.eqtl.gtex$pvalue)
    hist(-log10(tbl.eqtl.gtex$pvalue))
    table(tbl.eqtl.gtex$gene)
    if(!grepl("chr", tbl.eqtl.gtex$chrom[1]))
       tbl.eqtl.gtex$chrom <- paste0("chr", tbl.eqtl.gtex$chrom)

    title <- sprintf("gtex-eqtls-%s-%s-%d-%d-pval.ge02.RData",
                     targetGene, chrom, loc.start, loc.end)
    save(tbl.eqtl.gtex, file=title)

} # fetch.eqtls
#----------------------------------------------------------------------------------------------------
calculate.tf.target.correlations <- function()
{

  data.dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/expression"
  list.files(data.dir)
  filename <- "brandLabDifferentiationTimeCourse-27171x28-namesCorrected.RData"
  f <- file.path(data.dir, filename)
  file.exists(f)
  mtx.rna <- get(load(f))

  cors <- lapply(rownames(mtx.rna),
                 function(gene)
                    cor(mtx.rna[gene,], mtx.rna[targetGene,], method="spearman"))
  names(cors) <- rownames(mtx.rna)
  failures <- which(is.na(cors))
  length(failures)  # 3719
  cors <- cors[-failures]
  tbl.cor <- data.frame(gene=names(cors), cor=as.numeric(cors))
  feature.col.name <- "brandLab.rna.cor"

  ft$addGeneFeature(tbl.cor, feature.name=feature.col.name, default.value=0)

} # calculate.tf.target.correlations
#----------------------------------------------------------------------------------------------------
targetGene <- "ZBTB7A"
targetTissue <- "GTEx_V8.Whole_Blood"
tbl.fimo <- get(load("tbl.fimo.ZBTB7A.chr19:3255534-4745156.RData"))

ft <- FeatureTable$new(target.gene=targetGene, reference.genome="hg38")
ft$setFundamentalRegions(tbl.fimo)

tbl.eqtl <- get(load("gtex-eqtls-ZBTB7A-chr19-3987816-4132224-pval.ge02.RData"))
checkTrue(all(c("chrom", "rsid", "hg38", "pvalue") %in% colnames(tbl.eqtl)))
tbl.eqtl$start <- tbl.eqtl$hg38 - 1
tbl.eqtl$end <- tbl.eqtl$hg38

feature.guide <- list(blood.eqtl.rsid="rsid", blood.eqtl.pvalue="pvalue")
default.values <- list(blood.eqtl.rsid="", blood.eqtl.pvalue=1)
ft$addRegionFeature(tbl.eqtl, feature.genome="hg38", feature.guide, default.values)

calculate.tf.target.correlations()
tbl.ft <- ft$getTable()

tbl.ft.strong <- subset(tbl.ft, blood.eqtl.pvalue < 1e-8 &
                                abs(brandLab.rna.cor) > 0.6)
dim(tbl.ft.strong)
rsids.oi <- unique(tbl.ft.strong$blood.eqtl.rsid)
subset(tbl.eqtl, rsid %in% rsids.oi)

   #-------------------------------------------------------------
   #  save tbl.ft and tbl.eqtl for use in
   # ~/github/TrenaProjectErythropoiesis/explore/zbtb7a/zbtb7a.R
   #-------------------------------------------------------------

f <- file.path("~/github/TrenaProjectErythropoiesis/explore/zbtb7a",
               "tbl.ft.eqtl.RData")
save(tbl.ft, tbl.eqtl, file=f)
