library(MotifDb)
library(trena)
library(TrenaProjectHG38.generic)
library(igvR)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tp")){
   tp <- TrenaProjectHG38.generic()
   tbl.gene <- getTranscriptsTable(tp, "NFE2")
   tss <- tbl.gene$tss  # 54,295,760   this is the short transcript
   tss <- 54301037      # use this instead - the tss of the two longer transcripts
   start <- tss - 10000
   end   <- tss + 10000
   tbl.region <- data.frame(chrom="chr12", start=start, end=end)
   }

if(!exists("igv")){
   igv <- igvR()
   setGenome(igv, "hg38")
   showGenomicRegion(igv, "NFE2")
   tp <- TrenaProjectErythropoiesis()
   setTargetGene(tp, "NFE2")
   #tbl.gh <- getEnhancers(tp)
   #tbl.gh$width <- with(tbl.gh, 1 + end - start)
   #tbl.gh <- subset(tbl.gh, width < 5000)
   #track <- DataFrameQuantitativeTrack("gh", tbl.gh[, c(1,2,3,11)], color="blue", autoscale=FALSE, min=0, max=50)
   #displayTrack(igv, track)
   showGenomicRegion(igv, "chr12:54,289,590-54,311,542")
   showGenomicRegion(igv, "chr12:54,282,156-54,326,061")
   showGenomicRegion(igv, "chr12:54,290,497-54,315,164")
   with(tbl.region, showGenomicRegion(igv, sprintf("%s:%d-%d", chrom, start, end)))
   }
#------------------------------------------------------------------------------------------------------------------------
if(!exists("fimoBatch")){
   source("~/github/fimoService/batchMode/fimoBatchTools.R")   # works on hagfish & khaleesi
   meme.file <- "jaspar2018-hocomoco.meme"
   motifs <- query(MotifDb, "hsapiens", c("jaspar2018", "hocomoco"))
   length(motifs)  # 1177
   export(motifs, con=meme.file, format="meme")
   }
#------------------------------------------------------------------------------------------------------------------------
fimo.bioc.motifs <- function()
{
   tbl.fimoMotifs <- fimoBatch(tbl.region, 1e-3, "hg38", meme.file, quiet = TRUE)
   dim(tbl.fimoMotifs)
   tbl.fimoMotifs <- as.data.frame(gscores(phast.7, GRanges(tbl.fimoMotifs)), stringsAsFactors=FALSE)
   colnames(tbl.fimoMotifs)[grep("seqnames", colnames(tbl.fimoMotifs))] <- "chrom"
   colnames(tbl.fimoMotifs)[grep("default", colnames(tbl.fimoMotifs))] <- "phast7"

   mm <- MotifMatcher("hg38", as.list(motifs))
   tbl.biocMotifs <- findMatchesByChromosomalRegion(mm, tbl.region, 75)
   colnames(tbl.biocMotifs)[3:4] <- c("start", "end")

   dim(tbl.biocMotifs)
   tbl.biocMotifs.x <- as.data.frame(gscores(phast.7, GRanges(tbl.biocMotifs)), stringsAsFactors=FALSE)

   colnames(tbl.biocMotifs.x)[grep("seqnames", colnames(tbl.biocMotifs.x))] <- "chrom"
   colnames(tbl.biocMotifs.x)[grep("default", colnames(tbl.biocMotifs.x))] <- "phast7"
   tbl.biocMotifs <- tbl.biocMotifs.x
   tfs <- mcols(MotifDb[tbl.biocMotifs$motifName])$geneSymbol
   tbl.biocMotifs$tf <- tfs
   save(tbl.fimoMotifs, tbl.biocMotifs, file="tbls.motifs.RData")

} # bioc.motifs
#------------------------------------------------------------------------------------------------------------------------
select.tfs.by.motif.scores <- function()
{
   load("tbls.motifs.RData")
   fimo.score <- 1e-4
   phast.score <- 0.95
   bioc.score <- .90

   tfs.fimo <- subset(tbl.fimoMotifs, p.value < fimo.score & phast7 > phast.score)$tf
   tfs.bioc <- unique(subset(tbl.biocMotifs, motifRelativeScore > bioc.score & phast7 > phast.score)$tf)
   shared.tfs <- intersect(tfs.fimo, tfs.bioc)
   tfs.fimo.only <- setdiff(tfs.fimo, tfs.bioc)
   tfs.bioc.only <- setdiff(tfs.bioc, tfs.fimo)
   printf("shared: %6d   fimo.only: %6d  bioc.only: %6d",
          length(shared.tfs), length(tfs.fimo.only), length(tfs.bioc.only))

} # select.tfs.by.motif.scores
#------------------------------------------------------------------------------------------------------------------------
