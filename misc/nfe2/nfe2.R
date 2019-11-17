library(igvR)
library(TrenaProjectErythropoiesis)
library(GenomicScores)
library(phastCons7way.UCSC.hg38); phast.7 <- phastCons7way.UCSC.hg38
library(trenaSGM)
if(!exists("tbl.benchmark")){
   tbl.benchmark <- get(load(system.file(package="TrenaValidator", "extdata", "tbl.A.RData")))
   tbl.benchmark$pubmed.count <- unlist(lapply(strsplit(tbl.benchmark$pubmedID_from_curated_resources, ","), length))
   as.data.frame(t(subset(tbl.benchmark, TF=="GATA1" & target=="NFE2")))
   }

#  TF                                          GATA1
#  target                                       NFE2
#  effect                                          0
#  score                                           A
#  is_evidence_curateddatabase                  TRUE
#  is_evidence_chipSeq                         FALSE
#  is_evidence_TFbindingMotif                   TRUE
#  is_evidence_coexpression                     TRUE
#  which_curateddatabase               HTRIdb,trrust
#  which_chipSeq                                none
#  which_TFbindingMotif                 hocomoco_v11
#  which_coexpression                    ARACNe-GTEx
#  pubmedID_from_curated_resources 16648487,20564185
#  kegg_pathway                                    -
#  pubmed.count                                    2

#------------------------------------------------------------------------------------------------------------------------
# TF                                          GATA1
# target                                       NFE2
# effect                                          0
# score                                           A
# is_evidence_curateddatabase                  TRUE
# is_evidence_chipSeq                         FALSE
# is_evidence_TFbindingMotif                   TRUE
# is_evidence_coexpression                     TRUE
# which_curateddatabase               HTRIdb,trrust
# which_chipSeq                                none
# which_TFbindingMotif                 hocomoco_v11
# which_coexpression                    ARACNe-GTEx
# pubmedID_from_curated_resources 16648487,20564185
# kegg_pathway                                    -
# pubmed.count                                    2
#------------------------------------------------------------------------------------------------------------------------
# https://www.ncbi.nlm.nih.gov/pubmed/16648487
#
# Functional interaction of CP2 with GATA-1 in the regulation of erythroid promoters.   2006
# cp2 -> TFCP2
# Hsapiens-HOCOMOCOv10-TFCP2_HUMAN.H10MO.D
# Hsapiens-jaspar2018-TFCP2-MA0145.3

# We observed that binding sites for the ubiquitously expressed transcription factor CP2 [TFCP2] were
# present in regulatory regions of multiple erythroid genes. In these regions, the CP2 binding site
# was adjacent to a site for the erythroid factor GATA-1. Using three such regulatory regions (from
# genes encoding the transcription factors GATA-1, EKLF, and p45 NF-E2), we demonstrated the
# functional importance of the adjacent CP2/GATA-1 sites. In particular, CP2 binds to the GATA-1 HS2
# enhancer, generating a ternary complex with GATA-1 and DNA. Mutations in the CP2 consensus greatly
# impaired HS2 activity in transient transfection assays with K562 cells. Similar results were
# obtained by transfection of EKLF and p45 NF-E2 mutant constructs. Chromatin immunoprecipitation
# with K562 cells showed that CP2 binds in vivo to all three regulatory elements and that both
# GATA-1 and CP2 were present on the same GATA-1 and EKLF regulatory elements. Adjacent CP2/GATA-1
# sites may represent a novel module for erythroid expression of a number of genes. Additionally,
# coimmunoprecipitation and glutathione S-transferase pull-down experiments demonstrated a physical
# interaction between GATA-1 and CP2. This may contribute to the functional cooperation between
# these factors and provide an explanation for the important role of ubiquitous CP2 in the
# regulation of erythroid genes.
#
# 10778    chr12 54295889 54295904    16      + TFCP2 17.2360 1.56e-07 CCCTGCCTGGGCCAGA Hsapiens-HOCOMOCOv10-TFCP2_HUMAN.H10MO.D 0.95625
#
# lapply(query(MotifDb, c("sapiens", "TFCP2"), c("jaspar2018", "hocomoco")), consensusString)
# $`Hsapiens-jaspar2018-TFCP2-MA0145.3`
# [1] "AAACCGGTT?"
#
# $`Hsapiens-HOCOMOCOv10-TFCP2_HUMAN.H10MO.D`
# [1] "?CCTG??C??GCC?GA"   where the trailing GA is followed by TAAGA: perhaps the adjacency mentioned above
#
#------------------------------------------------------------------------------------------------------------------------
# https://www.ncbi.nlm.nih.gov/pubmed/20564185
#
# Over-expression of EDAG in the myeloid cell line 32D: induction of GATA-1 expression and erythroid/megakaryocytic
# phenotype.  (2010)
#
# Erythroid differentiation-associated gene (EDAG), a hematopoietic tissue-specific transcription
# regulator, plays a key role in maintaining the homeostasis of hematopoietic lineage
# commitment. However, the mechanism and genes regulated by EDAG remain unknown. In this study, we
# showed that overexpression of EDAG in a myeloid cell line 32D led to an erythroid phenotype with
# increased number of benzidine-positive cells and up-regulation of erythroid specific surface
# marker TER119. The megakaryocytic specific marker CD61 was also induced significantly. Using a
# genome-wide microarray analysis and a twofold change cutoff, we identified 332 genes with reduced
# expression and 288 genes with increased expression. Among up-regulation genes, transcription
# factor GATA-1 and its target genes including EKLF, NF-E2, Gfi-1b, hemogen, SCL, hemoglobin alpha,
# beta and megakaryocytic gene GPIX were increased. Silencing of EDAG by RNA interference in K562
# cells resulted in down-regulation of these genes. Taken together, EDAG functions as a positive
# regulator of erythroid/megakaryocytic differentiation in 32D cells associated with the induction
# of GATA-1 and its target genes.
#
#------------------------------------------------------------------------------------------------------------------------
# from google search: gata1 nf32
#
# Genetic Analysis of Hierarchical Regulation for Gata1 and NF-E2 p45 Gene Expression in Megakaryopoiesis
#  mouse, 2010. "a unique in vivo validation of the hierarchical relationship between GATA1 and p45 in megakaryocytes"
#
# Insight into GATA1 transcriptional activity through interrogation of cis elements disrupted in human erythroid disorders
# pnas (2016) very strong paper but seems only to discuss downstream effects of GATA1 and NFE2
#
# https://mcb.asm.org/content/25/4/1215
# minireview: GATA1 Function, a Paradigm for Transcription Factors in Hematopoiesis
# includes this suggestion that gata1 mutant causes significant decrease in nfe2 expression
# (p45 is an alias for nfe2.  nfe1 is an alias for gata1)

# (ii) MaFK and p45 NF-E2.In addition to the results seen with other megakaryocyte-specific genes,
# the expression of the transcription factors MafK and p45 NF-E2 is significantly decreased in
# megakaryocytes expressing an N-finger mutant of GATA1 (V205G) and in GATA1-deficient
# megakaryocytes (96, 116). p45 NF-E2 p45 and small Maf factors are critical for terminal
# differentiation of megakaryocytes (71, 100). This suggests that the attenuated expression of the
# essential transcription factors NF-E2 p45 and MafK is a major cause of the megakaryocytic
# phenotype of GATA1 mutations.


# Regulation of Mouse p45 NF-E2 Transcription by an Erythroid-specific GATA-dependent Intronic Alternative Promoter*
# mouse, 2000.   two alternate promoters.
#
# The erythroid-enriched transcription factor NF-E2 is composed of two subunits, p45 and p18, the
# former of which is mainly expressed in the hematopoietic system. We have isolated and
# characterized the mouse p45 NF-E2 gene; we show here that, similar to the human gene, the mouse
# gene has two alternative promoters, which are differentially active during development and in
# different hematopoietic cells. Transcripts from the distal promoter are present in both erythroid
# and myeloid cells; however, transcripts from an alternative proximal 1b promoter, lying in the
# first intron, are abundant in erythroid cells, but barely detectable in myeloid cells. During
# development, both transcripts are detectable in yolk sac, fetal liver, and bone
# marrow. Transfection experiments show that proximal promoter 1b has a strong activity in erythroid
# cells, which is completely dependent on the integrity of a palindromic GATA-1 binding site. In
# contrast, the distal promoter 1a is not active in this assay. When the promoter 1b is placed 3′ to
# the promoter 1a and reporter gene, in an arrangement that resembles the natural one, it acts as an
# enhancer to stimulate the activity of the upstream promoter la.

# previous two articles mentioned in plos 2011
#  Hierarchical Differentiation of Myeloid Progenitors Is Encoded in the Transcription Factor Network
#
# we examplarily discuss five such cases. (i) NF-E2 is regulated by GATA-1 and SCL (TAL1), but specifically
# important for megakaryocytic development [23]–[25].
#
if(!exists("igv")){
   igv <- igvR()
   setGenome(igv, "hg38")
   showGenomicRegion(igv, "NFE2")
   tp <- TrenaProjectErythropoiesis()
   setTargetGene(tp, "NFE2")
   tbl.gh <- getEnhancers(tp)
   tbl.gh$width <- with(tbl.gh, 1 + end - start)
   tbl.gh <- subset(tbl.gh, width < 5000)
   track <- DataFrameQuantitativeTrack("gh", tbl.gh[, c(1,2,3,11)], color="blue", autoscale=FALSE, min=0, max=50)
   displayTrack(igv, track)
   showGenomicRegion(igv, "chr12:54,289,590-54,311,542")
   showGenomicRegion(igv, "chr12:54,282,156-54,326,061")
   showGenomicRegion(igv, "chr12:54,290,497-54,315,164")
   }
#------------------------------------------------------------------------------------------------------------------------
round.numeric.columns.in.dataframe <- function(tbl, digits=2, pvalColumnNames="lassoPValue")
{
   tbl.pvals <- data.frame()
   tbl.main <- tbl

   if(!(all(is.na(pvalColumnNames)))){
     pval.cols <- grep(pvalColumnNames, colnames(tbl))
     stopifnot(length(pval.cols) == length(pvalColumnNames))
     tbl.pvals <- tbl[, pval.cols, drop=FALSE]
     tbl.main <- tbl[, -pval.cols, drop=FALSE]
     }
   numeric_columns <- sapply(tbl.main, mode) == 'numeric'
   tbl.main[numeric_columns] <-  round(tbl.main[numeric_columns], digits)
   if(ncol(tbl.pvals) > 0){
   tbl.pvals <- apply(tbl.pvals, 2, function(col) as.numeric(formatC(col, format = "e", digits = 2)))
   }

  tbl.out <- cbind(tbl.main, tbl.pvals)[, colnames(tbl)]
  tbl.out

   } # round.numeric.columns.in.dataframe

#------------------------------------------------------------------------------------------------------------------------
conservationTrack <- function()
{
   loc <- getGenomicRegion(igv)
   starts <- with(loc, seq(start, end, by=5))
   ends <- starts + 5
   count <- length(starts)
   tbl.blocks <- data.frame(chrom=rep(loc$chrom, count), start=starts, end=ends, stringsAsFactors=FALSE)
   tbl.cons7 <- as.data.frame(gscores(phast.7, GRanges(tbl.blocks)), stringsAsFactors=FALSE)
   tbl.cons7$chrom <- as.character(tbl.cons7$seqnames)
   tbl.cons7 <- tbl.cons7[, c("chrom", "start", "end", "default")]
   track <- DataFrameQuantitativeTrack("phast7", tbl.cons7, autoscale=TRUE, color="red")
   displayTrack(igv, track)

} # conservationTrack
#------------------------------------------------------------------------------------------------------------------------
fimoConservationTable <- function()
{
  source("~/github/fimoService/batchMode/fimoBatchTools.R")   # works on hagfish & khaleesi
  meme.file <- "jaspar2018-hocomoco.meme"
  motifs <- query(MotifDb, "hsapiens", c("jaspar2018", "hocomoco"))
  length(motifs)  # 1177
  export(motifs, con=meme.file, format="meme")

  roi <- getGenomicRegion(igv)
  tbl.regions <- with(roi, data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE))
  fimo.threshold <- 1e-5
  tbl.match <- fimoBatch(tbl.regions, matchThreshold=fimo.threshold, genomeName="hg38", pwmFile=meme.file)
  dim(tbl.match)

  tbl.matchCons <- as.data.frame(gscores(phast.7, GRanges(tbl.match)), stringsAsFactors=FALSE)
  dim(tbl.matchCons)
  tbl.matchCons <- subset(tbl.matchCons, default > 0.95)
  dim(tbl.matchCons)
  return(tbl.matchCons)

} # fimoConservationTable
#------------------------------------------------------------------------------------------------------------------------
demo_NFE2_models <- function()
{
   library(TrenaProjectLymphocyte)
   library(org.Hs.eg.db)
   tp <- TrenaProjectLymphocyte();

   genome <- "hg38"
   targetGene <- "NFE2"
   setTargetGene(tp, targetGene)
   tbl.info <- getTranscriptsTable(tp)

   chromosome <- tbl.info$chrom
   tss <- tbl.info$tss
      # strand-aware start and end: atf1 is on the + strand
   start <- tss - 50000
   end   <- tss + 50000

   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)
   file <- system.file(package="TrenaProjectLymphocyte", "extdata", "expression","GTEX.wholeBlood.rna-seq-geneSymbols.22330x407.RData")
   mtx.blood <- get(load(file))

   file <- system.file(package="TrenaProjectLymphocyte", "extdata", "expression",
                       "GTEX.lymphocyte.rna-seq-geneSymbols.21415x130.RData")
   mtx.lymphocyte <- get(load(file))

   file <- system.file(package="TrenaProjectErythropoiesis", "extdata", "expression", "brandLabDifferentiationTimeCourse-27171x28.RData")
   mtx.marjorie <- get(load(file))

   file <- "~/github/TrenaProjectErythropoiesis/prep/import/buenrostro/GSE74246_RNAseq_All_Counts.txt"
   file.exists(file)
   tbl <- read.table(file, sep="\t", as.is=TRUE, header=TRUE,nrow=-1)
   rownames(tbl) <- tbl[, 1]
   tbl <- tbl[, -1]
   mtx.buenrosto <- asinh(as.matrix(tbl))
   fivenum(mtx.buenrosto)



   tbl.matchCons <- fimoConservationTable()
   candidate.tfs <- unique(tbl.matchCons$tf)
   length(candidate.tfs)  # 57

   noDNA.recipe <- list(title="noDNA.matchCons",
                        type="noDNA.tfsSupplied",
                        matrix=mtx.marjorie,
                        #matrix=mtx.blood,
                        #matrix=mtx.lymphocyte,
                        #matrix=mtx.buenrosto,
                        candidateTFs=candidate.tfs,
                        tfPool=allKnownTFs(),
                        tfPrefilterCorrelation=0.2,
                        annotationDbFile=dbfile(org.Hs.eg.db),
                        orderModelByColumn="rfScore",
                        solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman", "xgboost"),
                        quiet=TRUE)
   builder <- NoDnaModelBuilder(genome, targetGene,  noDNA.recipe, quiet=TRUE)
   x <- build(builder)
   tbl.model <- x$model[order(abs(x$model$pearsonCoeff), decreasing=TRUE),]
   tbl.tfbs.counts <- as.data.frame(sort(table(tbl.matchCons$tf)))
   bindingSiteCount <- merge(tbl.model, tbl.tfbs.counts, by.x="gene", by.y="Var1")$Freq
   tbl.model$bindingSites <- bindingSiteCount

   tbl.model.strong <- subset(tbl.model, abs(pearsonCoeff) > 0.5)
   displayBindingSites(tbl.model.strong, tbl.matchCons)
   mtx.model <- as.matrix(tbl.model.strong[, -1])
   rownames(mtx.model) <- tbl.model.strong$gene

   tbl.model.trimmed <- as.data.frame(round.numeric.columns.in.dataframe(mtx.model))
   save(tbl.model.trimmed, file="brand.tbl.model.trimmed")



   recipe.all <- noDNA.recipe
   recipe.all$candidateTFs <- allKnownTFs()
   builder <- NoDnaModelBuilder(genome, targetGene,  recipe.all, quiet=TRUE)
   x1 <- build(builder)



      #----------------------------------------------------------------------------------------------------
      # first, build a model with "placenta2", an early version of the placenta footprint database
      #----------------------------------------------------------------------------------------------------

   recipe <- list(title="NFE2",
                  type="footprint.database",
                  regions=tbl.regions,
                  geneSymbol=targetGene,
                  tss=tss,
                  matrix=mtx.blood,
                  db.host="khaleesi.systemsbiology.net",
                  db.port=5432,
                  databases=list("lymphoblast_hint_16", "lymphoblast_hint_20"),
                  annotationDbFile=dbfile(org.Hs.eg.db),
                  motifDiscovery="builtinFimo",
                  tfPool=allKnownTFs(),
                  tfMapping="MotifDB",
                  tfPrefilterCorrelation=0.1,
                  orderModelByColumn="rfScore",
                  solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman", "xgboost"))

   fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene,  recipe, quiet=FALSE)
   x <- build(fpBuilder)

   recipe$matrix <- mtx.marjorie
   fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene,  recipe, quiet=FALSE)
   x2 <- build(fpBuilder)


   file <- "~/github/TrenaProjectErythropoiesis/prep/import/buenrostro/GSE74246_RNAseq_All_Counts.txt"
   file.exists(file)
   tbl <- read.table(file, sep="\t", as.is=TRUE, header=TRUE,nrow=-1)
   rownames(tbl) <- tbl[, 1]
   tbl <- tbl[, -1]
   mtx <- asinh(as.matrix(tbl))
   fivenum(mtx)

   suppressWarnings(
      logTimingInfo("placenta2 db, +/- 5bk generic promoter on ATF1", system.time(x2 <- build(fpBuilder)))
      )

} # demo_NFE2_models
#------------------------------------------------------------------------------------------------------------------------
displayBindingSites <- function(tbl.model, tbl.matchCons)
{
   tfs <- subset(tbl.model, abs(pearsonCoeff) > 0.5)$gene

   for(one.tf in tfs){
      tbl.bs <- subset(tbl.matchCons, tf==one.tf)[, c("seqnames", "start", "end", "matched_sequence")]
      colnames(tbl.bs) <- c("chrom", "start", "end", "seq")
      tbl.bs$chrom <- as.character(tbl.bs$chrom)
      track <- DataFrameAnnotationTrack(one.tf, tbl.bs, color="random", trackHeight=25, displayMode="EXPANDED")
      displayTrack(igv, track)
      }

} # displayBindingSites
#------------------------------------------------------------------------------------------------------------------------
getATACseq <- function()
{
   roi <- getGenomicRegion(igv)
   chromosome <- roi$chrom
   start.loc <- roi$start
   end.loc <- roi$end

   directory <- "~/github/TrenaProjectErythropoiesis/prep/import/atacPeaks"
   files <- grep("narrowPeak$", list.files(directory), value=TRUE)
   result <- list()

   for(file in files){
      full.path <- file.path(directory, file)
      track.name <- sub("_hg38_macs2_.*$", "", sub("ATAC_Cord_", "", file))
      tbl.atac <- read.table(full.path, sep="\t", as.is=TRUE)
      colnames(tbl.atac) <- c("chrom", "start", "end", "name", "c5", "strand", "c7", "c8", "c9", "c10")
      tbl.atac.region <- subset(tbl.atac, chrom==chromosome & start >= start.loc & end <= end.loc)
      if(nrow(tbl.atac.region) > 0){
         tbl.atac.region$sample <- track.name
         result[[track.name]] <- tbl.atac.region
         }
      } # files

   tbl.out <- do.call(rbind, result)
   rownames(tbl.out) <- NULL

   tbl.out

} # getATACseq
#------------------------------------------------------------------------------------------------------------------------
displayATACseq <- function()
{
   library (RColorBrewer)
   totalColorCount <- 12
   colors <- brewer.pal(8, "Dark2")
   currentColorNumber <- 0

   tbl.all <- getATACseq()
   samples <- unique(tbl.all$sample)
   current.day.string <- ""
   color <- colors[1]

   for(current.sample in samples){
      this.day.string <- strsplit(current.sample, "_")[[1]][1]
      if(this.day.string != current.day.string){
         currentColorNumber <- (currentColorNumber %% totalColorCount) + 1
         color <- colors[currentColorNumber]
         current.day.string <- this.day.string
         }
      tbl.atac.sub <- subset(tbl.all, sample == current.sample)
      track.name <- current.sample
      track <- DataFrameQuantitativeTrack(track.name, tbl.atac.sub[, c("chrom", "start", "end", "c10")],
                                          color, autoscale=FALSE, min=0, max=430, trackHeight=30)
      displayTrack(igv, track)
      } # for samples

   tbl.regions.condensed <- as.data.frame(union(GRanges(tbl.all[, c("chrom", "start", "end")]),
                                                GRanges(tbl.all[, c("chrom", "start", "end")])))[, c("seqnames", "start", "end")]
   colnames(tbl.regions.condensed) <- c("chrom", "start", "end")
   tbl.regions.condensed$chrom <- as.character(tbl.regions.condensed$chrom)
   lapply(tbl.regions.condensed, class)

   #state$tbl.regions.condensed <- tbl.regions.condensed
   track <- DataFrameAnnotationTrack("atac combined", tbl.regions.condensed, color="black")
   displayTrack(igv, track)

} # displayATACseq
#------------------------------------------------------------------------------------------------------------------------
# oddly, no methylation data in the immediate vicinity of nfe2.
addMethylationTracks <- function()
{
  library(AnnotationHub)
  ah <- AnnotationHub()
  ah.human <- subset(ah, species == "Homo sapiens")
  histone.tracks <- query(ah.human, c("H3K4me3", "Gm12878", "Peak", "narrow"))  # 3 tracks
  descriptions <- histone.tracks$description
  titles <- histone.tracks$title
  colors <- rep(terrain.colors(6), 4)
  color.index <- 0

  tbl.roi <- as.data.frame(getGenomicRegion(igv), stringsAsFactors=FALSE)

  for(i in seq_len(length(histone.tracks))){
     name <- names(histone.tracks)[i]
     color.index <- color.index + 1
     gr <- histone.tracks[[name]]
     ov <- findOverlaps(gr, GRanges(tbl.roi))
     roi.histones <- gr[queryHits(ov)]
     track.histones <- GRangesQuantitativeTrack(titles[i], roi.histones[, "pValue"],
                                              color=colors[color.index], trackHeight=50,
                                              autoscale=TRUE)
     displayTrack(igv, track.histones)
     } # for track

} # addMethylationTracks
#------------------------------------------------------------------------------------------------------------------------
