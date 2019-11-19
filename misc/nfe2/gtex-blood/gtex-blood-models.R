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
build.with.all.GO.tfs <- function()
{
   tfs.all <- sort(unique(select(org.Hs.eg.db, keys="GO:0003700", keytype="GOALL", columns="SYMBOL")$SYMBOL))
   length(tfs.all)  # 1663

   target.gene <- "NFE2"
   tfs <- intersect(tfs.all, rownames(mtx.blood.lps))
   length(tfs)  # 1652

   solver <- EnsembleSolver(mtx.blood.lps, target.gene, tfs, geneCutoff=1.0)
   tbl <- run(solver)
   dim(tbl)
   new.order <- order(abs(tbl$pearsonCoeff), decreasing=TRUE)
   tbl <- tbl[new.order,]
   rownames(tbl) <- NULL
   tbl.lps.goAll <- tbl # [1:100,]
   match(c("GATA1", "TAL1", "KLF1"), tbl.lps.goAll$gene)   # 371 725 824



} # build.with.all.GO.tfs
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
   tbl.lps.motifs <- tbl # [1:100,]

   match(c("GATA1", "TAL1", "KLF1"), tbl.lps.motifs$gene)   # 125, 225, 252



   motifdb.tfs <- mcols(query(MotifDb, c("sapiens"), c("jaspar2018", "hocomoco")))$geneSymbol
   tfs <- intersect(motifdb.tfs, rownames(mtx.blood.lps))
   length(tfs)
      # no.motifs <- setdiff(head(tbl.lps$gene, n=20), motifdb.tfs)
      # length(no.motifs)
      # match(no.motifs, tbl.lps$gene)

   solver <- EnsembleSolver(mtx.blood.lps, target.gene, tfs, geneCutoff=1.0)
   tbl.lps <- run(solver)
   dim(tbl.lps)
   new.order <- order(abs(tbl.lps$pearsonCoeff), decreasing=TRUE)
   tbl.lps <- tbl.lps[new.order,]
   match(c("GATA1", "TAL1", "KLF1"), tbl.lps$gene)   # 125 225 252

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
build.with.fimo.and.phast <- function()
{
   load("../tbls.motifs.RData")
   fimo.score <- 1e-5
   phast.score <- 0.90

   tbl.tfs.elite <- subset(tbl.fimoMotifs, p.value <= fimo.score & phast7 >= phast.score)
   dim(tbl.tfs.elite)
   tfs.elite <- sort(unique(tbl.tfs.elite$tf))
   length(tfs.elite)  # 52
   match(c("GATA1", "TAL1", "KLF1"), tfs.elite)   # 13 41 16


   solver <- EnsembleSolver(mtx.blood.lps, target.gene, tfs.elite, geneCutoff=1.0, solverNames=solverNames)
   tbl <- run(solver)
   dim(tbl)
   new.order <- order(abs(tbl$pearsonCoeff), decreasing=TRUE)
   tbl <- tbl[new.order,]
   rownames(tbl) <- NULL
   tbl.fimo.phast.stringent <- tbl
   match(c("GATA1", "TAL1", "KLF1"), tbl.fimo.phast.stringent$gene)   # 15 27 30


   tbl.tfs.weak <- subset(tbl.fimoMotifs, p.value <= 1e-4 & phast7 > 0.2)
   dim(tbl.tfs.weak)   # 2846
   tfs.weak <- unique(tbl.tfs.weak$tf)
   length(tfs.weak)  # 525
   match(c("GATA1", "TAL1", "KLF1"), tfs.weak)   # 281 282 330

   solver <- EnsembleSolver(mtx.blood.lps, target.gene, tfs.weak, geneCutoff=1.0, solverNames=solverNames)
   tbl <- run(solver)
   dim(tbl)  # 379 7
   new.order <- order(abs(tbl$pearsonCoeff), decreasing=TRUE)
   tbl <- tbl[new.order,]
   rownames(tbl) <- NULL
   tbl.fimo.phast.weak <- tbl
   match(c("GATA1", "TAL1", "KLF1"), tbl.fimo.phast.weak$gene)   #  95 168 190

   tbl.tfbs.counts <- as.data.frame(sort(table(tbl.fimo.strong$tf)))
   bindingSiteCount <- merge(tbl.blood.lps.fimo, tbl.tfbs.counts, by.x="gene", by.y="Var1")$Freq
   tbl.blood.lps.fimo$tfbs.strong <- bindingSiteCount

   tbl.tfbs.counts <- as.data.frame(sort(table(tbl.fimo.weak$tf)))
   bindingSiteCount <- merge(tbl.blood.lps.fimo, tbl.tfbs.counts, by.x="gene", by.y="Var1")$Freq
   tbl.blood.lps.fimo$tfbs.weak <- bindingSiteCount

   match(c("GATA1", "TAL1", "KLF1"), tbl.blood.lps.fimo$gene)   # 8 4 5

   save(tbl.blood.lps.fimo, file="trena.model.fimo.phast7.RData")




} # build.with.fimo.and.phast
#------------------------------------------------------------------------------------------------------------------------
