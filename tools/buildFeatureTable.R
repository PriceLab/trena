library(trena)
library(ghdb)
library(RUnit)
library(RPostgreSQL)
library(EndophenotypeExplorer)

library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifbreakR)
library(MotifDb)
library(BiocParallel)

#----------------------------------------------------------------------------------------------------
buildFT = R6Class("buildFT",

    #--------------------------------------------------------------------------------
    private = list(targetGene=NULL,
                   gtex.tissue=NULL,
                   ft=NULL,
                   rosmap.eqtls=NULL,
                   fimo.file=NULL,
                   chrom=NULL,
                   start=NULL,
                   end=NULL,
                   etx=NULL,
                   mtx.rna=NULL,
                   tbl.filtered=NULL,
                   tbl.filtered.collapsed=NULL,
                   tbl.trena=NULL,
                   tbl.breaks=NULL,
                   tbl.compositeScores=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(

        initialize = function(targetGene, gtex.tissue, fimo.file,
                              study.region.start=NA, study.region.end=NA){
           private$targetGene <- targetGene
           private$gtex.tissue <- gtex.tissue
           private$ft <- FeatureTable$new(target.gene=targetGene, reference.genome="hg39")
           private$etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD")
           stopifnot(file.exists(fimo.file))
           private$fimo.file <- fimo.file
           tbl.fimo <- get(load(fimo.file))
           if(!is.na(study.region.start))
               tbl.fimo <- subset(tbl.fimo, start >= study.region.start)
           if(!is.na(study.region.end))
               tbl.fimo <- subset(tbl.fimo, end <= study.region.end)
           private$start <- min(tbl.fimo$start)
           private$end   <- max(tbl.fimo$end)
           private$chrom <- tbl.fimo$chrom[1]
           p.value.col <- grep("^p\\.value$", colnames(tbl.fimo))
           if(length(p.value.col) == 1)
               colnames(tbl.fimo)[p.value.col] <- "fimo.pval"
           score.col <- grep("^score$", colnames(tbl.fimo))
           if(length(score.col) == 1)
               colnames(tbl.fimo)[score.col] <- "fimo.score"
           private$ft$setFundamentalRegions(tbl.fimo)
           },

       #------------------------------------------------------------
       getGtexTissue = function(){
           return(private$ft$gtex.tissue)
           },

       #------------------------------------------------------------
       getTable = function(){
           return(private$ft$getTable())
           },

       #------------------------------------------------------------
       getRnaMatrixCodes = function(){
           names(private$etx$get.rna.matrix.codes())
           },

       #------------------------------------------------------------
       getRnaMatrix = function(){
           invisible(private$mtx.rna)
           },

       #------------------------------------------------------------
       setFilteredTable = function(tbl.filtered){
           private$tbl.filtered <- tbl.filtered
           },

       #------------------------------------------------------------
       getFilteredTable = function(){
           invisible(private$tbl.filtered)
           },

       #------------------------------------------------------------
       getCollapsedFilteredTable = function(){
           invisible(private$tbl.filtered.collapsed)
           },

       #------------------------------------------------------------
       getTrenaTable = function(){
           private$tbl.trena
           },

       #------------------------------------------------------------
       # a lengthy computation. this function can sidestep it
       setBreaksTable = function(tbl.breaks){
           private$tbl.breaks <- tbl.breaks
           },

       #------------------------------------------------------------
       getBreaksTable = function(){
           private$tbl.breaks
           },

       #------------------------------------------------------------
       getRegion = function(){
           return(list(chrom=private$chrom,
                       start=private$start,
                       end=private$end,
                       width=1 + private$end - private$start))
           },

       #------------------------------------------------------------
       add.rosmap.eqtls = function(){
           if(file.exists("rosmap.eqtls.RData")){
               printf("adding rosmap eqtls from prior run")
               tbl.eqtl <- get(load("rosmap.eqtls.RData"))
           } else {
              printf("adding rosmap eqtls from genereg2021 on khaleesi")
              geneRegDB <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="genereg2021", host="khaleesi")
              query <- sprintf("select * from eqtl2 where chrom='%s' and hg38 > %d and hg38 < %d  and genesymbol ='%s';",
                               private$chrom, private$start, private$end, private$targetGene)
              tbl.eqtl <- dbGetQuery(geneRegDB, query)
              tbl.eqtl$start <- tbl.eqtl$hg38 - 1
              tbl.eqtl$end <- tbl.eqtl$hg38
              }
           if(!"score" %in% colnames(tbl.eqtl))
               tbl.eqtl$score <- with(tbl.eqtl, -log10(pvalue) * beta)
           feature.guide <- list(rosmap.eqtl.rsid="rsid", rosmap.eqtl.pvalue="pvalue",
                                 rosmap.eqtl.beta="beta", rosmap.eqtl.score="score")
           checkTrue(all(as.character(feature.guide) %in% colnames(tbl.eqtl)))
           default.values <- list(rosmap.eqtl.rsid="", rosmap.eqtl.pvalue=1,
                                  rosmap.eqtl.beta=0, rosmap.eqtl.score=0)
           private$ft$addRegionFeature(tbl.eqtl, feature.genome="hg38", feature.guide, default.values)
           private$rosmap.eqtls <- tbl.eqtl
           }, # add.rosmap.eqtls

       #------------------------------------------------------------
       add.gtex.eqtls = function(gtex.eqtl.file){
           stopifnot(file.exists(gtex.eqtl.file))
           tbl.eqtl <- get(load(gtex.eqtl.file))
           tbl.eqtl.sub <- subset(tbl.eqtl, hg38 >= private$start & hg38 <= private$end &
                                            chrom==private$chrom & id == private$gtex.tissue &
                                            gene==private$targetGene)
           tbl.eqtl.sub$start <- tbl.eqtl.sub$hg38 - 1
           tbl.eqtl.sub$end <- tbl.eqtl.sub$hg38
           if(!"score" %in% colnames(tbl.eqtl.sub))
               tbl.eqtl.sub$score <- with(tbl.eqtl.sub, -log10(pvalue) * beta)
           feature.guide <- list("rsid", "pvalue", "beta", "score")
           feature.names <- sprintf("gtex.eqtl.%s", feature.guide)
           names(feature.guide) <- feature.names
           default.values <- list("", 1, 0, 0)
           names(default.values) <- feature.names
           private$ft$addRegionFeature(tbl.eqtl.sub, feature.genome="hg38", feature.guide, default.values)
           }, # add.gtex.eqtls

       #------------------------------------------------------------
       add.mayo.atac = function(){

           data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
           file <- "mayoAllPeaks.merged.96064x4.RData"
           full.path <- file.path(data.dir, file)
           stopifnot(file.exists(full.path))
           tbl.atac <- get(load(full.path))
           tbl.atac.sub <- subset(tbl.atac, start >= private$start & end <= private$end &
                                            chrom==private$chrom)
           tbl.atac.sub$status <- TRUE
           feature.guide <-  list(mayoAtac="status")
           default.values <- list(mayoAtac=FALSE)
           private$ft$addRegionFeature(tbl.atac.sub, feature.genome="hg38", feature.guide, default.values)
           }, # add.mayo.atac

       #------------------------------------------------------------
       add.boca.atac = function(){

           data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
           file <- "boca-hg38-consensus-ATAC.RData"
           full.path <- file.path(data.dir, file)
           stopifnot(file.exists(full.path))
           tbl.atac <- get(load(full.path))
           tbl.atac.sub <- subset(tbl.atac, start >= private$start & end <= private$end &
                                            chrom==private$chrom)
           tbl.atac.sub$status <- TRUE
           feature.guide <-  list(bocaAtac="status")
           default.values <- list(bocaAtac=FALSE)
           private$ft$addRegionFeature(tbl.atac.sub, feature.genome="hg38", feature.guide, default.values)
           }, # add.boca.atac

       #------------------------------------------------------------
       add.expression.correlation = function(matrix.code, short.name){
           mtx <- private$etx$get.rna.matrix(matrix.code)
           private$mtx.rna <- mtx
           cor.func <- function(gene){
              cor(mtx[private$targetGene,], mtx[gene, ], method="spearman", use="pairwise.complete")
              }
           x <- lapply(rownames(mtx), cor.func)
           tbl.cor <- data.frame(gene=rownames(mtx), cor=as.numeric(x))
           new.colname <- sprintf("rna.cor.%s", short.name)
           colnames(tbl.cor)[2] <- new.colname
           private$ft$addGeneFeature(tbl.cor, feature.name=new.colname, default.value=0)
           },

       #------------------------------------------------------------
       add.genehancer = function(tissue){

           stopifnot(tissue %in% c("brain", "all"))
           ghdb <- GeneHancerDB()
           tbl.gh <- retrieveEnhancersFromDatabase(ghdb, private$targetGene, tissues=tissue)
           feature.name <- sprintf("gh.%s", tissue)
           feature.guide <- list("combinedscore")
           default.values <- list(0)
           names(feature.guide) <- feature.name
           names(default.values) <- feature.name
           private$ft$addRegionFeature(tbl.gh, feature.genome="hg38", feature.guide, default.values)
           }, # add.genehancer

       #------------------------------------------------------------
       run.trena = function(tfs, add.rsids.column=FALSE){
          solver <- EnsembleSolver(private$etx$get.rna.matrix(private$gtex.tissue),
                                   targetGene=private$targetGene,
                                   candidateRegulators=tfs,
                                   solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest"), #, "xgboost"),
                                   geneCutoff=1.0)
          tbl.trena <- run(solver)
          new.order <- order(tbl.trena$rfScore, decreasing=TRUE)
          tbl.trena <- tbl.trena[new.order,]
          tbl.fft <-  self$getCollapsedFilteredTable()
          tbl.trena$tfbs <- as.integer(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.fft, tf==gene))))
          gtex.eqtl.column.name <- "gtex.eqtl.rsid"
          x <- lapply(tbl.trena$gene, function(gene) subset(tbl.fft, tf==gene)[, ..gtex.eqtl.column.name])
          rsids <- as.character(lapply(x, function(rsids) paste(unlist(rsids), collapse=";")))
          if(add.rsids.column)
             tbl.trena$rsids <- rsids
          rownames(tbl.trena) <- NULL
          tbl.trena$betaLasso <- round(tbl.trena$betaLasso, digits=3)
          tbl.trena$betaRidge <- round(tbl.trena$betaRidge, digits=3)
          tbl.trena$spearmanCoeff <- round(tbl.trena$spearmanCoeff, digits=3)
          tbl.trena$pearsonCoeff <- round(tbl.trena$pearsonCoeff, digits=3)
          tbl.trena$rfScore <- round(tbl.trena$rfScore, digits=3)
          private$tbl.trena <- tbl.trena
          invisible(private$tbl.trena)
          }, # run.trena

       #------------------------------------------------------------
          # if a tf has overlapping binding sites, and the same rsid, reduce it to a single tfbs
       reduce.tbl.filtered = function(quiet=TRUE){
          stopifnot(!is.null(private$tbl.filtered))
          tbl.fft <- private$tbl.filtered
          tf.xtab <- as.list(table(tbl.fft$tf))
          tfs.mult <- names(tf.xtab[tf.xtab > 1])
          tfs.single <- names(tf.xtab[tf.xtab == 1])
          tbl.fftC <- subset(tbl.fft, tf %in% tfs.single)  # "C":  collapsed
          gtex.eqtl.column.name <- "gtex.eqtl.rsid"
          for(tf.mult in tfs.mult){
              tbl.sub <- subset(tbl.fft, tf==tf.mult)
              all.rsids <- as.character(unlist(tbl.sub[, ..gtex.eqtl.column.name]))
                 # match returns only the first hit
              keeper.rows <- match(unique(all.rsids), all.rsids)
              tbl.sub.collapsed <- tbl.sub[keeper.rows,]
              if(!quiet) printf("adding %d rows for %s", nrow(tbl.sub.collapsed), tf.mult)
              tbl.fftC <- rbind(tbl.fftC, tbl.sub.collapsed)
              }
          #browser()
          private$tbl.filtered.collapsed <- tbl.fftC
          xyz <- 99
          },
       #------------------------------------------------------------
       findMotifBreaks = function(pval.threshold){
              # we only care about breaks in tfbs motifs in the final
          tbl.fftC <- self$getCollapsedFilteredTable()
          stopifnot("rosmap.eqtl.pvalue" %in% colnames(tbl.fftC))
          tbl.fftC.sub <- subset(tbl.fftC, rosmap.eqtl.pvalue <= pval.threshold & nchar(rosmap.eqtl.rsid) > 0)
          if(nrow(tbl.fftC.sub) == 0)
              return(data.frame())
          rsids.sub <- unique(tbl.fftC.sub$rosmap.eqtl.rsid)
          motifs.oi <- unique(tbl.fftC.sub$motif_id)
          printf("----- snps.from.rsid, starting.  snp count: %d  motif count: %d", length(rsids.sub), length(motifs.oi))
          print(load("../shared/1883.rosmap.snps.from.fimo.region.for.motifbreakR.RData"))
          snps.gr.sub <- subset(snps.gr, SNP_id %in% rsids.sub)
          #times <- system.time(snps.gr <- snps.from.rsid(rsid=rsids.sub,
          #                                               dbSNP=SNPlocs.Hsapiens.dbSNP155.GRCh38,
          #                                               search.genome=BSgenome.Hsapiens.UCSC.hg38))
          #printf("----- snps.from.rsid, elapsed: %6.3f", times[["elapsed"]])
          motifs <- MotifDb[motifs.oi]

          printf("----- motifbreakR, starting")
          times <- system.time(
          results <- motifbreakR(snpList = snps.gr.sub,
                                 filterp = TRUE,
                                 pwmList = motifs,
                                 show.neutral=FALSE,
                                 method = c("ic", "log", "notrans")[1],
                                 bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                 BPPARAM = BiocParallel::bpparam(),
                                 verbose=TRUE))
          printf("----- motifbreakR, elapsed: %6.3f", times[["elapsed"]])
          tbl.breaks <- as.data.frame(results, row.names=NULL)
          if(nrow(tbl.breaks) == 0) return()
          tbl.breaks <- subset(tbl.breaks, effect=="strong")
          if(nrow(tbl.breaks) == 0) return()
          tbl.breaks$pctDelta <- round(with(tbl.breaks, pctAlt-pctRef), digits=2)
          new.order <- order(abs(tbl.breaks$pctDelta), decreasing=TRUE)
          tbl.breaks <- tbl.breaks[new.order,]
          dups <- which(duplicated(tbl.breaks[, c("SNP_id", "geneSymbol")]))
          if(length(dups) > 0)
              tbl.breaks <- tbl.breaks[-dups,]
          rownames(tbl.breaks) <- NULL
          private$tbl.breaks <- tbl.breaks
          }, # findMotifBreaks

       #------------------------------------------------------------
       calculateCompositeScores = function(debug=FALSE, quiet=TRUE){

          tbl.trena.oi <- head(private$tbl.trena, n=-1)
          tbl.trena.oi$rfNorm <- tbl.trena.oi$rfScore/max(tbl.trena.oi$rfScore)
          tbl.trena.oi$betaLassoNorm <- abs(tbl.trena.oi$betaLasso)/max(abs(tbl.trena.oi$betaLasso))
          tbl.fftC <- private$tbl.filtered.collapsed

          tbl.rosmap.snps <- subset(tbl.fftC, abs(rosmap.eqtl.score) > 0 & tf %in% tbl.trena.oi$gene)
          tbl.breaks <- self$getBreaksTable()
          dups <- which(duplicated(tbl.breaks[, c("SNP_id", "geneSymbol")]))
          if(length(dups) > 0)
              tbl.breaks <- tbl.breaks[-dups,]
          dim(tbl.breaks)

            #------------------------------------------------------------
            # now score the impact of these breaks, the sum of breakage
            # scores for each tf:
            # motifBreak.score <- with(tbl.tf, abs(pctDelta) * 100)
            #    eqtl.score <- with(tbl.tf, -log10(gtex.eqtl.pval)* abs(gtex.eqtl.beta) * 100)
            #  trena.score <- with(tbl.tf, (abs(betaLasso) * 100) + (rfNorm * 10))
            #   tfbs.score <- 1/tfbs.count
            # score <- trena.score * eqtl.score * motifBreak.score * tfbs.score
            #
            #   abs(pctDelta) * tf's (rfNorm  + abs(betaLasso))
            #
            #-----------------------------------------------------------

           composite.scores <- list()

           for(TF in tbl.trena.oi$gene){
               tf.breakage.score <- 0
               tbl.breaks.sub <- subset(tbl.breaks, geneSymbol==TF & pctDelta < 0)
               for(rsid in tbl.breaks.sub$SNP_id){
                   tbl.trena.sub <- subset(tbl.trena.oi, gene==TF)
                   motifbreakR.pctDelta <- abs(subset(tbl.breaks.sub, SNP_id==rsid)$pctDelta) * 10
                   trena.rfNorm <- tbl.trena.sub$rfNorm * 10
                   trena.betaLasso.score <- tbl.trena.sub$betaLassoNorm * 10
                   trena.score <- trena.betaLasso.score + trena.rfNorm
                   rosmap.eqtl.score <- abs(subset(tbl.rosmap.snps, rosmap.eqtl.rsid==rsid & tf==TF)$rosmap.eqtl.score * 2)
                   if(length(rosmap.eqtl.score) == 0)
                       rosmap.eqtl.score <- 0
                   increment <- rosmap.eqtl.score * motifbreakR.pctDelta * trena.score
                   if(!quiet){
                       printf("--- %s (%s)", TF, rsid)
                       printf("   pctDelta: %5.2f", motifbreakR.pctDelta)
                       printf("     rfNorm: %5.2f", trena.rfNorm)
                       printf("  betaLasso: %5.2f", trena.betaLasso.score)
                       printf(" trenaScore: %5.2f", trena.score)
                       printf("  rosmap eqtl: %5.2f", rosmap.eqtl.score)
                       printf("  increment: %5.2f", increment)
                       }
                   tf.breakage.score <- tf.breakage.score + increment
                   } # for rsid
               composite.score <- as.integer(tf.breakage.score)
               composite.scores[[TF]] <- composite.score
               printf(" ====== %8s: %6d", TF, composite.score)
               } # for tf
          private$tbl.trena$composite.score <- 0
          indices <- match(names(composite.scores), tbl.trena.oi$gene)
          tbl.trena.oi$composite.score[indices] <- unlist(composite.scores, use.names=FALSE)
          private$tbl.trena <- tbl.trena.oi
          } # calculateCompositeScores

       ) # public


    ) # class
#----------------------------------------------------------------------------------------------------
