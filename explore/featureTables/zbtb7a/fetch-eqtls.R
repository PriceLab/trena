targetGene <- "ZBTB7A"

   require(ebi.eqtls)
   ee <- ebi.eqtls$new()
   tbl.cat <- ee$getCatalog()
   toi <- subset(tbl.cat, Tissue_group %in% c("blood", "Whole Blood") &
                          quant_method=="ge")$unique_id
      #  "Lepik_2017.blood" "TwinsUK.blood"    "GTEx.blood" "GTEx_V8.Whole_Blood"

   start.small <- 4058249
   end.small <-   4062002
   start.144k <- 3987816
   end.144k <-   4132224
   start.long <- 2652636
   end.long   <- 4719616

   chrom <- "chr19"
   loc.start <- start.144k
   loc.end   <- end.144k

   loc.start <- start.long
   loc.end <- end.long

   loc.start <- start.small
   loc.end <- end.small

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

