library(trena)
library(ghdb)


ghdb <- GeneHancerDB()
targetGene <- "NDUFS2"
tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
tbl.gh$score <- asinh(tbl.gh$combinedscore)

data.dir <- "data"

   #-----------------------------------------------------------------------------------
   # create the FeatureTable, initialize with fimo results trimmed to GeneHancer range
   #----------------------------------------------------------------------------------

ft <- FeatureTable$new(target.gene=targetGene, reference.genome="hg38")
tbl.fimo <- get(load(file.path(data.dir, "tbl.fimo.NDUFS2.RData")))
dim(tbl.fimo)   # 2117704       9

                # trim to quasi-informative limits provided by GeneHancer
tbl.fimo <- subset(tbl.fimo, start > (min(tbl.gh$start) - 1000) &
                             end   < (max(tbl.gh$end)   + 1000))
dim(tbl.fimo)   #  828415      9
tbl.track <- tbl.fimo[, c("chrom", "start", "end", "p.value")]
tbl.track$p.value <- -log10(tbl.track$p.value)

ft$setFundamentalRegions(tbl.fimo)

   #-------------------------------------------------------------------
   # add GTEx brain eQTLs for NDUFS2, cortex only
   #-------------------------------------------------------------------

tbl.eqtl.gtex <- get(load(file.path(data.dir, "gtex.brain.eqtls.chr1-160684126-161688825.RData")))
dim(tbl.eqtl.gtex)   # 1690734     8
tbl.eqtl.gtex <- subset(tbl.eqtl.gtex, gene==targetGene)
dim(tbl.eqtl.gtex)   #   29553     8
table(tbl.eqtl.gtex$id)
tbl.eqtl.gtex <- subset(tbl.eqtl.gtex, id=="GTEx_V8.Brain_Cortex" & pvalue < 0.10)
dim(tbl.eqtl.gtex)   #     391     9
tbl.eqtl.gtex$score <- -log10(tbl.eqtl.gtex$pvalue)
tbl.eqtl.gtex$start <- tbl.eqtl.gtex$hg38 - 1
tbl.eqtl.gtex$end <- tbl.eqtl.gtex$hg38

tbl.eqtl.gtex.ready <- tbl.eqtl.gtex[, c("chrom", "start", "end", "pvalue", "score", "rsid")]

feature.guide <- list(gtex.eqtl.rsid="rsid", gtex.eqtl.pvalue="pvalue")
default.values <- list(gtex.eqtl.rsid="", gtex.eqtl.pvalue=1)
ft$addRegionFeature(tbl.eqtl.gtex.ready, feature.genome="hg38", feature.guide, default.values)


   #-------------------------------------------------------------------
   # add rosmap brain eQTLs for NDUFS2
   #-------------------------------------------------------------------

tbl.eqtl.rosmap <- get(load(file.path(data.dir, "rosmap.eqtls.RData")))
dim(tbl.eqtl.rosmap)   # 1894 14
tbl.eqtl.rosmap <- subset(tbl.eqtl.rosmap, pvalue < 0.10)
dim(tbl.eqtl.rosmap)   #  346 14
tbl.eqtl.rosmap.ready <- tbl.eqtl.rosmap[, c("chrom", "hg38", "hg38", "pvalue", "rsid")]
tbl.eqtl.rosmap.ready$score <- -log10(tbl.eqtl.rosmap$pvalue)
tbl.eqtl.rosmap.ready$start <- tbl.eqtl.rosmap$hg38 - 1
tbl.eqtl.rosmap.ready$end <- tbl.eqtl.rosmap$hg38

feature.guide <- list(rosmap.eqtl.rsid="rsid", rosmap.eqtl.pvalue="pvalue")
default.values <- list(rosmap.eqtl.rsid="", rosmap.eqtl.pvalue=1)
ft$addRegionFeature(tbl.eqtl.rosmap.ready, feature.genome="hg38", feature.guide, default.values)

   #-------------------------------------------------------------------
   # add mayo atac
   #-------------------------------------------------------------------

filename <- "mayoAllPeaks.merged.96064x4.RData"
tbl.mayoAtac <- get(load(file.path(data.dir, filename)))
dim(tbl.mayoAtac)
tbl.mayoAtac <- subset(tbl.mayoAtac, chrom==tbl.eqtl.rosmap$chrom[1] &
                                     start >= min(tbl.eqtl.rosmap$hg38) - 1000 &
                                     end   <= max(tbl.eqtl.rosmap$hg38) + 1000)
dim(tbl.mayoAtac)  # 36 5

tbl.mayoAtac$status <- TRUE   # we don't have scores, alas
feature.guide <- list(mayo.atac="status")
default.values <- list(mayo.atac=FALSE)
ft$addRegionFeature(tbl.mayoAtac, feature.genome="hg38", feature.guide, default.values)
# tbl.ft <- ft$getTable()

   #-------------------------------------------------------------------
   # add the small set of eQTL european haplotype scores from haploreg
   #-------------------------------------------------------------------

# f <- system.file(package="trena", "extdata", "featureTable", "tbl.haploreg.european.rs4575098.RData")
tbl.haplo <- read.table(file.path(data.dir, "haploreg-rs4575098-0.2.tsv"), sep="\t", header=TRUE)
dim(tbl.haplo)
tbl.haplo$start <- tbl.haplo$hg38 - 1
tbl.haplo$end   <- tbl.haplo$hg38

feature.guide <-  list(haploreg.rsid="rsid", haploreg.rSquared="rSquared")
default.values <- list(haploreg.rsid="",     haploreg.rSquared=0)
ft$addRegionFeature(tbl.haplo, feature.genome="hg38", feature.guide, default.values)

   #-------------------------------------------------------------------------
   # add NDUFS2/TF correlated expression scores from GTEx Brain_Cortex.RData
   #-------------------------------------------------------------------------

mtx <- get(load(file.path(data.dir, "Brain_Cortex.RData")))
dim(mtx)   # 24821 206


genes.to.consider <- c(targetGene, unique(intersect(tbl.fimo$tf, rownames(mtx))))
length(genes.to.consider)
gene.correlations <- lapply(genes.to.consider, function(gene) cor(mtx[targetGene,], mtx[gene,], method="spearman", use="complete"))
tbl.cor <- data.frame(gene=genes.to.consider, cor=as.numeric(gene.correlations))
dim(tbl.cor)  # 354 2

subset(tbl.cor, abs(cor) > 0.66)
   #       gene        cor
   # 1   NDUFS2  1.0000000
   # 113 TCF7L1 -0.6641626
   # 162   RXRB  0.7452202
   # 283   NFIA -0.6654400
   # 314  FOXK2  0.6602555

feature.guide <- list(gtex.cortex.rna.cor="cor")
default.values <- list(gtex.cortex.rna.cor=0)
ft$addGeneFeature(tbl.cor, feature.name="gtex.cortex.rna.cor", default.value=0)
tbl.ft <- ft$getTable()


   #-------------------------------------------------------------------------
   # add simulated enhancer promoter interactions
   #-------------------------------------------------------------------------

# library(GenomicRanges)
# library(annotatr)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# annotation.sought <- "hg38_genes_promoters"
# annotations <- build_annotations(genome="hg38", annotations=annotation.sought)
# gr.atac <- GRanges(tbl.mayoAtac)
# tbl.ov <- as.data.frame(findOverlaps(gr.atac, annotations))
#
#
# tbl.atac.promoters <- tbl.mayoAtac[tbl.ov[,1],]
# keepers <- tbl.ov[,2]
# tbl.atac.promoters$gene <- annotations[keepers]$symbol
# rownames(tbl.atac.promoters) <- NULL
# tbl.atac.promoters <- unique(subset(tbl.atac.promoters, gene==targetGene))
# dim(tbl.atac.promoters) #   1 7 -> of all the mayo atac, only  chr1:161201843-16120272
#                         #  is annotated as an NDUFS2 promoter
# with(tbl.atac.promoters, 1 + end - start)  # 878


# generate some random enhancer/promoter links which
#    originate in mayo atac open chromatin
#    terminate within this 879 bp promoter
# how large are these regions reported to be?
#
# tbl.pe <- read.table("~/github/igvR/inst/extdata/gm12878_loops.bedpe.gz", sep="\t")
# colnames(tbl.pe) <- c("chr1", "start1", "end1", "chr1", "start2", "end2")
# with(tbl.pe, fivenum(end1-start1))  # 5000  5000  5000 10000 10000
# with(tbl.pe, fivenum(end2-start2))  # 5000  5000  5000 10000 10000
#
#   # 	Targeted high-resolution chromosome conformation capture at genome-wide scale in mouse erythroid cells
# tbl.view <- read.table("GSE160229_Capture_Viewpoints.txt.gz", sep="\t", nrow=-1)
# dim(tbl.view) # 7197    4
# head(tbl.view)
# colnames(tbl.view) <- c("gene", "chrom", "start", "end")
# tbl.view$end <- as.numeric(tbl.view$end)
#
# fivenum(tbl.view$end-tbl.view$start)  # 69.0       461.0       745.0      1167.5 207101398.0
# order(table(tbl.view$end-tbl.view$start))
#    # from this, a leap:  imaginge that actual binding sites as 70 bp
#
#
# tbl.atac.promoters
#      chrom     start       end width strand status   gene
#       chr1 161201843 161202721   879      *   TRUE NDUFS2
# create 20 landing sites in the promoter

# starts <- sample(161201843:(161202721-70), size=nrow(tbl.mayoAtac))
# tbl.landing <- data.frame(chrom="chr1", start=starts, end=(starts+70))
# starts <- unlist(lapply(seq_len(nrow(tbl.mayoAtac)),
#        function(r) {
#            start <- tbl.mayoAtac[r, "start"];
#            end   <- tbl.mayoAtac[r,"end"]
#            sample(start: end-70, size=1)
#            }))
# tbl.launch <- data.frame(chrom="chr1", start=starts, endd=starts+70)
# tbl.pe <- cbind(tbl.launch, tbl.landing)
# colnames(tbl.pe) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
# dim(tbl.pe)
# random.subset <- sample(seq_len(nrow(tbl.pe)), size=5)
# tbl.pe <- tbl.pe[random.subset,]
# rownames(tbl.pe) <- NULL
#
tbl.pe <- get(load(file.path(data.dir, "tbl.pairedEnd.small.NDUFS2.RData")))

    #------------------------------------------------------------
    # any tfbs which is on either end of these paired end relations
    # deserves to be annotated to reflect that.
    # first the "launch" sites (enhancers) then the "landing" sites
    # (promoters)
    #------------------------------------------------------------

feature.guide <- list(launch.3c="status")
default.values <- list(launch.3c=FALSE)
tbl.launch <- tbl.pe[, c("chrom1", "start1", "end1")]
colnames(tbl.launch) <- c("chrom", "start", "end")
tbl.launch$status <- TRUE
ft$addRegionFeature(tbl.launch, feature.genome="hg38", feature.guide, default.values)

feature.guide <- list(landing.3c="status")
default.values <- list(landing.3c=FALSE)
tbl.landing <- tbl.pe[, c("chrom2", "start2", "end2")]
colnames(tbl.landing) <- c("chrom", "start", "end")
tbl.landing$status <- TRUE
ft$addRegionFeature(tbl.landing, feature.genome="hg38", feature.guide, default.values)


tbl.ft <- ft$getTable()
dim(tbl.ft)   # 2117704      19
tfs <- unique(tbl.ft$tf)
length(tfs)  # 535

  # tbl.fft: FilteredFeatureTable

tbl.fft <- subset(tbl.ft, gtex.eqtl.pvalue < 0.01 & rosmap.eqtl.pvalue < 0.01)
tfs <- unique(tbl.fft$tf)
length(tfs)  # 266


solvers <- c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest")

solver <- EnsembleSolver(mtx,  targetGene, tfs, solverNames=solvers)
system.time(tbl.trena <- run(solver))

tbl.trena$tfbs <- as.integer(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.fft, tf==gene))))
gtex.eqtl.column.name <- "gtex.eqtl.rsid"

if(any(grepl(gtex.eqtl.column.name, colnames(tbl.fft)))){
   rsids.raw  <- lapply(tbl.trena$gene, function(gene) subset(tbl.fft, tf==gene)[, ..gtex.eqtl.column.name])
   rsids.by.tf <- as.character(lapply(rsids.raw, function(rsids) paste(unlist(rsids), collapse=";")))
   tbl.trena$rsids <- rsids.by.tf
   }

new.order <- order(tbl.trena$rfScore, decreasing=TRUE)
tbl.trena <- tbl.trena[new.order,]
rownames(tbl.trena) <- NULL
tbl.trena$betaLasso <- round(tbl.trena$betaLasso, digits=3)
tbl.trena$betaRidge <- round(tbl.trena$betaRidge, digits=3)
tbl.trena$spearmanCoeff <- round(tbl.trena$spearmanCoeff, digits=3)
tbl.trena$pearsonCoeff <- round(tbl.trena$pearsonCoeff, digits=3)
tbl.trena$rfScore <- round(tbl.trena$rfScore, digits=3)

options(width=1000)
print(tbl.trena)
