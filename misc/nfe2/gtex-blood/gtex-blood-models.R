library(trena)
library(rnaSeqNormalizer)
library(GTEx)
library(trenaSGM)
library(MotifDb)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("gtex"))
    gtex <- GTEx()
if(!exists("mtx.blood.counts")){
   mtx.blood.counts <- createSubMatrix(gtex, "blood")
   normalizer.lps <- rnaSeqNormalizer.gtex(mtx.blood.counts, algorithm="log+scale", duplicate.selection.statistic="median")
   mtx.blood.lps <- getNormalizedMatrix(normalizer.lps)
   normalizer <- rnaSeqNormalizer.gtex(mtx.blood.counts, algorithm="asinh", duplicate.selection.statistic="median")
   mtx.blood.asinh <- getNormalizedMatrix(normalizer)
   }

#------------------------------------------------------------------------------------------------------------------------
build.with.tfs.with.motifs <- function()
{
      # use the log+scale normalizer

   fivenum(mtx.blood.lps)
   target.gene <- "NFE2"
   tfs.with.motifs <- sort(unique(mcols(query(MotifDb, c("sapiens"), c("jaspar2018", "hocomoco")))$geneSymbol))
   length(tfs.with.motifs)  # 780
   tfs <- intersect(tfs.with.motifs, rownames(mtx.blood.lps))
   length(tfs)  # 556

   solver <- EnsembleSolver(mtx.blood.lps, target.gene, tfs, geneCutoff=1.0)
   tbl <- run(solver)
   dim(tbl)
   new.order <- order(abs(tbl$pearsonCoeff), decreasing=TRUE)
   tbl <- tbl[new.order,]
   rownames(tbl) <- NULL
   tbl.lps <- tbl # [1:100,]

   match(c("GATA1", "TAL1", "KLF1"), tbl.lps$gene)   # 125, 225, 252

   motifdb.tfs <- mcols(query(MotifDb, c("sapiens"), c("jaspar2018", "hocomoo")))$geneSymbol
   no.motifs <- setdiff(head(tbl.lps$gene, n=20), motifdb.tfs)
   length(no.motifs)
   match(no.motifs, tbl.lps$gene)

     # use the asinh normalizer

   fivenum(mtx.blood.asinh)

   solver <- EnsembleSolver(mtx.blood.asinh, target.gene, tfs)
   tbl.asinh <- run(solver)
   new.order <- order(abs(tbl.asinh$pearsonCoeff), decreasing=TRUE)
   tbl.asinh <- tbl.asinh[new.order,]
   match(c("GATA1", "TAL1", "KLF1"), tbl.asinh$gene)   # NA 33 34

   save(tbl.lps, tbl.asinh, file="trena.models.noDNA")

} # build.with.tfs.with.motifs
#------------------------------------------------------------------------------------------------------------------------
build.with.tfs.with.matching.fimo.motifs <- function()
{
      # use the log+scale normalizer

   fivenum(mtx.blood.lps)
   target.gene <- "NFE2"

   load("../tbls.motifs.RData")
   fimo.score <- 1e-4
   phast.score <- 0.95

   tbl.fimo.strong <- subset(tbl.fimoMotifs, p.value <= fimo.score & phast7 >= phast.score)
   dim(tbl.fimo.strong)
   tfs <- sort(unique(tbl.fimo.strong$tf))
   length(tfs)  # 209

   tbl.fimo.weak <- subset(tbl.fimoMotifs, p.value <= 1e-4 & phast7 >= 0.5)
   dim(tbl.fimo.weak)

   solver <- EnsembleSolver(mtx.blood.lps, target.gene, tfs, geneCutoff=1.0)
   tbl <- run(solver)
   length(tbl$gene)
   dim(tbl)
   new.order <- order(abs(tbl$pearsonCoeff), decreasing=TRUE)
   tbl <- tbl[new.order,]
   rownames(tbl) <- NULL
   tbl.lps.fimo <- tbl
   dim(tbl.lps.fimo)

   tfs.in.model <- tbl.lps.fimo$gene
   length(tfs.in.model)
   tbl.fimo.strong.filtered <- subset(tbl.fimo.strong, tf %in% tfs.in.model)
   length(unique(tbl.fimo.strong.filtered$tf))

   tbl.tfbs.counts <- as.data.frame(sort(table(tbl.fimo.strong.filtered$tf)))
   dim(tbl.tfbs.counts)

   bindingSiteCount <- merge(tbl.lps.fimo, tbl.tfbs.counts, by.x="gene", by.y="Var1")$Freq
   length(bindingSiteCount)
   dim(tbl.lps.fimo)
   tbl.lps.fimo$tfbs.strong <- bindingSiteCount

   tbl.fimo.weak.filtered <- subset(tbl.fimo.weak, tf %in% tfs.in.model)
   tbl.tfbs.counts <- as.data.frame(sort(table(tbl.fimo.weak.filtered$tf)))
   bindingSiteCount <- merge(tbl.lps.fimo, tbl.tfbs.counts, by.x="gene", by.y="Var1", all.x=TRUE)$Freq
   tbl.lps.fimo$tfbs.weak <- bindingSiteCount

   match(c("GATA1", "TAL1", "KLF1"), tbl.lps.fimo$gene)   # 43 88 94

   save(tbl.lps.fimo, file="trena.models.fimo.phast7.RData")


} # build.with.tfs.with.matching.fimo.motifs
#------------------------------------------------------------------------------------------------------------------------

