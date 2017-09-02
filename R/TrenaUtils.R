.TrenaUtils <- setClass ("TrenaUtils",
                        representation = representation(
                           genomeName="character")
                        )
#------------------------------------------------------------------------------------------------------------------------
setGeneric('getRegulatoryChromosomalRegions',  signature='obj',
           function(obj, chromosome, chromStart, chromEnd, regulatoryRegionSources, targetGene, targetGeneTSS,
                    combine=FALSE, quiet=FALSE) standardGeneric("getRegulatoryChromosomalRegions"))


setGeneric('getRegulatoryTableColumnNames',  signature='obj', function(obj) standardGeneric ('getRegulatoryTableColumnNames'))
setGeneric('getGeneModelTableColumnNames',  signature='obj', function(obj) standardGeneric ('getGeneModelTableColumnNames'))
setGeneric('createGeneModel', signature='obj', function(obj, targetGene,  solverNames, tbl.regulatoryRegions, mtx)
              standardGeneric('createGeneModel'))
setGeneric('buildMultiModelGraph', signature='obj', function(obj, models) standardGeneric('buildMultiModelGraph'))
#setGeneric('expandRegulatoryRegionsTableByTF', signature='obj', function(obj, tbl.reg) standardGeneric('expandRegulatoryRegionsTableByTF'))
#setGeneric('addGeneModelLayout', signature='obj', function(obj, g, xPos.span=1500) standardGeneric('addGeneModelLayout'))
setGeneric('assessSnp', signature='obj', function(obj, pfms, variant, shoulder, pwmMatchMinimumAsPercentage, genomeName="hg38")
              standardGeneric('assessSnp'))
#------------------------------------------------------------------------------------------------------------------------
# a temporary hack: some constants
genome.db.uri <- "postgres://bddsrds.globusgenomics.org/hg38"   # has gtf and motifsgenes tables
#------------------------------------------------------------------------------------------------------------------------
TrenaUtils = function(genomeName, quiet=TRUE)
{
   stopifnot(genomeName %in% c("hg19", "hg38", "mm10"))

   obj <- .TrenaUtils(genomeName=genomeName)

   obj

} # constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod('getRegulatoryTableColumnNames', 'TrenaUtils',

      function(obj){
         c("chrom", "motifStart", "motifEnd", "motifName", "strand", "score", "length", "distance.from.tss", "id", "tf")
         })

#------------------------------------------------------------------------------------------------------------------------
setMethod('getGeneModelTableColumnNames', 'TrenaUtils',

      function(obj){
         c("tf", "randomForest", "pearson", "spearman", "betaLasso", "pcaMax", "concordance")
         })

#------------------------------------------------------------------------------------------------------------------------
.callFootprintFilter <- function(obj, source, chromosome, chromStart, chromEnd, targetGene, targetGeneTSS)
{
    chromLocString <- sprintf("%s:%d-%d", chromosome, chromStart, chromEnd)
    fpFilter <- FootprintFilter(genome.db.uri, source,  geneCenteredSpec=list(), regionsSpec=chromLocString)
    x.fp <- getCandidates(fpFilter)
    tbl.fp <- x.fp$tbl
    if(nrow(tbl.fp) == 0){
       warn("no footprints found in %s:%d-%d, targetGene is %s", chromosome, chromeStart, chromEnd, targetGene);
       return(tbl.fp)
       }

    colnames(tbl.fp) <- c("chrom", "motifStart", "motifEnd", "motifName", "length", "strand", "score1", "score", "score3", "tf")

    distance <- tbl.fp$motifStart - targetGeneTSS
    direction <- rep("upstream", length(distance))
    direction[which(distance < 0)] <- "downstream"
    tbl.fp$distance.from.tss <- distance
    tbl.fp$id <- sprintf("%s.fp.%s.%06d.%s", targetGene, direction, abs(distance), tbl.fp$motifName)
    tbl.fp <- tbl.fp[, getRegulatoryTableColumnNames(obj)]

    tbl.fp

} # .callFootprintFilter
#------------------------------------------------------------------------------------------------------------------------
.callHumanDHSFilter <- function(obj, chromosome, chromStart, chromEnd, targetGene, targetGeneTSS)
{
    printf("--- in .callHumanDHS")
    chromLocString <- sprintf("%s:%d-%d", chromosome, chromStart, chromEnd)
    dhsFilter <- HumanDHSFilter(genome="hg38",
                                encodeTableName="wgEncodeRegDnaseClustered",
                                pwmMatchPercentageThreshold=85L,
                                geneInfoDatabase.uri=genome.db.uri,
                                geneCenteredSpec=list(),
                                regionsSpec=chromLocString)

    x.dhs <- getCandidates(dhsFilter)

    if(all(is.na(x.dhs))){
       return(data.frame())
       }

    tbl.dhs <- x.dhs$tbl
    tbl.dhs$length <- nchar(tbl.dhs$match)
    distance <- tbl.dhs$motifStart - targetGeneTSS
    direction <- rep("upstream", length(distance))
    direction[which(distance < 0)] <- "downstream"

    colnames(tbl.dhs)[grep("motifRelativeScore", colnames(tbl.dhs))] <- "score"
    colnames(tbl.dhs)[grep("tfs", colnames(tbl.dhs))] <- "tf"
    tbl.dhs$distance.from.tss <- distance
    tbl.dhs$id <- sprintf("%s.dhs.%s.%06d.%s", targetGene, direction, abs(distance), tbl.dhs$motifName)

    tbl.dhs <- tbl.dhs[, getRegulatoryTableColumnNames(obj)]

    tbl.dhs

} # .callHumanDHSFilter
#------------------------------------------------------------------------------------------------------------------------
setMethod('getRegulatoryChromosomalRegions', 'TrenaUtils',

    function(obj, chromosome, chromStart, chromEnd, regulatoryRegionSources, targetGene, targetGeneTSS,
             combine=FALSE, quiet=FALSE){

         tbl.combined <- data.frame()
         result <- list()
           # some bookeeeping to permit duplicate sources, useful only in testing
         source.count <- 0
         all.source.names <- regulatoryRegionSources

         encodeDHS.source.index <- grep("encodeHumanDHS", regulatoryRegionSources)

         if(length(encodeDHS.source.index)){
            source.count <- source.count + 1
            regulatoryRegionSources <- regulatoryRegionSources[-encodeDHS.source.index]
            if(!quiet) printf("about to callHumanDHSFilter");
            tbl.dhs <- .callHumanDHSFilter(obj, chromosome, chromStart, chromEnd, targetGene, targetGeneTSS)
            result[[source.count]] <- tbl.dhs
            if(combine)
               tbl.combined <- rbind(tbl.combined, tbl.dhs)
            } # if encode DSH source requested


         for(source in regulatoryRegionSources){
            source.count <- source.count + 1
            if(!quiet) printf("about to call footprintFilter with source = '%s'", source);
            tbl.fp <- .callFootprintFilter(obj, source, chromosome, chromStart, chromEnd, targetGene, targetGeneTSS)
            if(combine)
               tbl.combined <- rbind(tbl.combined, tbl.fp)
            result[[source.count]] <- tbl.fp
            } # for source
         names(result) <- all.source.names
         if(combine)
            result[["all"]] <- tbl.combined
         result
         }) # getRegulatoryChromosomalRegions

#------------------------------------------------------------------------------------------------------------------------
#setMethod('expandRegulatoryRegionsTableByTF', 'TrenaUtils',
#
#     function(obj, tbl.reg){
#        tbl.trimmed <- subset(tbl.reg, nchar(tf) != 0)
#        tfs.split <- strsplit(tbl.trimmed$tf, ";")
#        #length(tfs.split) # [1] 36929
#        counts <- unlist(lapply(tfs.split, length))
#        tfs.split.vec <- unlist(tfs.split)
#        tbl.expanded <- expandRows(tbl.trimmed, counts, count.is.col=FALSE, drop=FALSE)
#        stopifnot(length(tfs.split.vec) == nrow(tbl.expanded))
#        tbl.expanded$tf <- tfs.split.vec
#        tbl.expanded
#        }) # expandRegulatoryRegionsTableByTF
#
#------------------------------------------------------------------------------------------------------------------------
setMethod('createGeneModel', 'TrenaUtils',

      function(obj, targetGene, solverNames, tbl.regulatoryRegions, mtx){

         #stopifnot(solvers=="randomForest")  # more solvers to come
         tfs <- sort(unique(unlist(strsplit(tbl.regulatoryRegions$tf, ";"))))
         tfs <- intersect(tfs, rownames(mtx))
         printf("tf candidate count: %d", length(tfs))

         #solver.wt <- RandomForestSolver(mtx, targetGene=targetGene, candidateRegulators=tfs)
         solver <- EnsembleSolver(mtx, targetGene=targetGene, candidateRegulators=tfs, solverNames)
         tbl.model  <-run(solver)
         #count <- nrow(model.wt$edges)
         #tbl.model <- data.frame(tf=rownames(model.wt$edges),
         #                        randomForest=model.wt$edges$IncNodePurity,
         #                        pearson=model.wt$edges$gene.cor,
         #                        spearman=rep(0, count),
         #                        betaLasso=rep(0, count),
         #                        pcaMax=rep(0, count),
         #                        concordance=rep(0, count),
         #                        stringsAsFactors=FALSE)
         tbl.model
      }) # createGeneModel

#------------------------------------------------------------------------------------------------------------------------
setMethod('buildMultiModelGraph', 'TrenaUtils',

  function (obj, models){

    g <- graphNEL(edgemode = "directed")
    model.names <- names(models)

    node.attribute.specs <- list(type="undefined",
                                 label="default node label",
                                 distance=0,
                                 pearson=0,
                                 randomForest=0,
                                 pcaMax=0,
                                 concordance=0,
                                 betaLasso=0,
                                 motif="",
                                 xPos=0,
                                 yPos=0)
    edge.attribute.spec <- list(edgeType="undefined")
    attribute.classes <- c("", model.names)  # "" (no prefix) is the currently displayed set of attibutes

      # create current version of these attributes, and then
      # per-model versions, which get mapped to current
      # in response to user's interactive choice on the cyjs user interface
      # the "current version" is, e.g., "distance".
      # per-model ("wt" and "mut" versions) become "wt.distance" and "mut.distance"
      # and are used by copying e.g. all wt.xxx attributes into the current (non-prefixed)
      # attribute, upon which the cyjs style is defined

    for(class.name in attribute.classes){
       class.name.prefix <- class.name  # with possible "." appended, permits standard and model-specific attributes
       if(nchar(class.name) > 0)
          class.name.prefix <- sprintf("%s.", class.name)
       noa.names.without.prefix <- names(node.attribute.specs)
       noa.names <- sprintf("%s%s", class.name.prefix, noa.names.without.prefix)
       noa.count <- length(node.attribute.specs)
       for(i in 1:noa.count){
          nodeDataDefaults(g, attr=noa.names[i]) <- node.attribute.specs[[noa.names.without.prefix[i]]]
          }
       } # for class

    edgeDataDefaults(g, attr = "edgeType") <- "undefined"

    tfs <- c()
    regulatoryRegions <- c()

    for(model in models){  # collect all the tf and regulatory region nodes
       tbl.model <- model$tbl.geneModel
       tfs <- unique(c(tfs, tbl.model$tf))
       tbl.reg <- model$tbl.regulatoryRegions
       regulatoryRegions <- unique(c(regulatoryRegions, tbl.reg$id))
       } # for model

    all.nodes <- unique(c(obj@targetGene, tfs, regulatoryRegions))
    g <- addNode(all.nodes, g)

    nodeData(g, obj@targetGene, "type") <- "targetGene"
    nodeData(g, tfs, "type")         <- "TF"
    nodeData(g, regulatoryRegions, "type")  <- "regulatoryRegion"
    nodeData(g, all.nodes, "label")  <- all.nodes

      # add edges, edge attribute, and the constant attributes for all of the regulatoryRegion nodes

    for(model in models){
       tfs <- model$tbl.regulatoryRegions$tf
       regRegions <- model$tbl.regulatoryRegions$id
       suppressWarnings(g <- addEdge(tfs, regRegions, g))
       edgeData(g,  tfs, regRegions, "edgeType") <- "bindsTo"
       suppressWarnings(g <- addEdge(regRegions, obj@targetGene, g))
       edgeData(g, regRegions, obj@targetGene, "edgeType") <- "regulatorySiteFor"
       nodeData(g, tbl.reg$id, "label") <- tbl.reg$motifName
       nodeData(g, tbl.reg$id, "distance") <- tbl.reg$distance.from.tss
       nodeData(g, tbl.reg$id, "motif") <- tbl.reg$motifName
       } # for model

      # now copy in the first model's tf node data

    tbl.model <- models[[1]]$tbl.geneModel
    nodeData(g, tbl.model$tf, attr="randomForest") <- tbl.model$randomForest
    nodeData(g, tbl.model$tf, attr="pearson") <- tbl.model$pearson

     # now copy in each of the model's tf node data in turn
    model.names <- names(models)
    for(model.name in model.names){
       tbl.model <- models[[model.name]]$tbl.geneModel
       noa.name <- sprintf("%s.%s", model.name, "randomForest")
       nodeData(g,  tbl.model$tf, attr=noa.name) <- tbl.model$randomForest
       noa.name <- sprintf("%s.%s", model.name, "pearson")
       nodeData(g,  tbl.model$tf, attr=noa.name) <- tbl.model$pearson
      } # for model.name

    g

    }) # buildMultiModelGraph

#------------------------------------------------------------------------------------------------------------------------
#setMethod('addGeneModelLayout', 'TrenaUtils',
#
#  function (obj, g, xPos.span=1500){
#    printf("--- addGeneModelLayout")
#    all.distances <- sort(unique(unlist(nodeData(g, attr='distance'), use.names=FALSE)))
#    print(all.distances)
#
#    fp.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "regulatoryRegion")]
#    tf.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "TF")]
#    targetGene.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "targetGene")]
#
#     # add in a zero in case all of the footprints are up or downstream of the 0 coordinate, the TSS
#    span.endpoints <- range(c(0, as.numeric(nodeData(g, fp.nodes, attr="distance"))))
#    span <- max(span.endpoints) - min(span.endpoints)
#    footprintLayoutFactor <- 1
#    printf("initial:  span: %d  footprintLayoutFactor: %f", span, footprintLayoutFactor)
#
#    footprintLayoutFactor <- xPos.span/span
#
#    #if(span < 600)  #
#    #   footprintLayoutFactor <- 600/span
#    #if(span > 1000)
#    #   footprintLayoutFactor <- span/1000
#
#    printf("corrected:  span: %d  footprintLayoutFactor: %f", span, footprintLayoutFactor)
#
#    xPos <- as.numeric(nodeData(g, fp.nodes, attr="distance")) * footprintLayoutFactor
#    yPos <- 0
#    nodeData(g, fp.nodes, "xPos") <- xPos
#    nodeData(g, fp.nodes, "yPos") <- yPos
#
#    adjusted.span.endpoints <- range(c(0, as.numeric(nodeData(g, fp.nodes, attr="xPos"))))
#    printf("raw span of footprints: %d   footprintLayoutFactor: %f  new span: %8.0f",
#           span, footprintLayoutFactor, abs(max(adjusted.span.endpoints) - min(adjusted.span.endpoints)))
#
#    tfs <- names(which(nodeData(g, attr="type") == "TF"))
#
#    for(tf in tfs){
#       footprint.neighbors <- edges(g)[[tf]]
#       if(length(footprint.neighbors) > 0){
#          footprint.positions <- as.integer(nodeData(g, footprint.neighbors, attr="xPos"))
#          new.xPos <- mean(footprint.positions)
#          if(is.na(new.xPos)) browser()
#          if(is.nan(new.xPos)) browser()
#          #printf("%8s: %5d", tf, new.xPos)
#          }
#       else{
#          new.xPos <- 0
#          }
#       nodeData(g, tf, "yPos") <- sample(300:1200, 1)
#       nodeData(g, tf, "xPos") <- new.xPos
#       } # for tf
#
#    nodeData(g, targetGene.nodes, "xPos") <- 0
#    nodeData(g, targetGene.nodes, "yPos") <- -200
#
#    g
#
#    }) # addGeneModelLayout
#
#------------------------------------------------------------------------------------------------------------------------
setMethod('assessSnp', 'TrenaUtils',

     function(obj, pfms=list(), variant, shoulder, pwmMatchMinimumAsPercentage, genomeName="hg38"){

        motifMatcher <- MotifMatcher(name=variant, genomeName=genomeName, pfms=pfms, quiet=TRUE)
        tbl.variant <- trena:::.parseVariantString(motifMatcher, variant)
        tbl.regions <- data.frame(chrom=tbl.variant$chrom,
                                  start=tbl.variant$loc-shoulder,
                                  end=tbl.variant$loc+shoulder,
                                  stringsAsFactors=FALSE)
        x.wt  <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions,
                                                pwmMatchMinimumAsPercentage=pwmMatchMinimumAsPercentage)
        if(nrow(x.wt$tbl) == 0){
           warning(sprintf("no motifs found in reference sequence in neighborhood of %s with shoulder %d",
                           variant, shoulder))
           return(data.frame())
           }


        x.mut <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions,
                                                pwmMatchMinimumAsPercentage=pwmMatchMinimumAsPercentage,
                                                variant=variant)



        if(nrow(x.mut$tbl) == 0){
           warning(sprintf("no motifs altered by %s with shoulder %d", variant, shoulder))
           return(data.frame())
           }

        tbl.wt.50 <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, 50)$tbl
        tbl.wt.50$signature <- sprintf("%s;%s;%s", tbl.wt.50$motifName, tbl.wt.50$motifStart, tbl.wt.50$strand)
        tbl.mut.50 <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, 50, variant=variant)$tbl
        tbl.mut.50$signature <- sprintf("%s;%s;%s", tbl.mut.50$motifName, tbl.mut.50$motifStart, tbl.mut.50$strand)

        tbl <- rbind(x.wt$tbl[, c(1,12, 2,3,4,5,7,8,13)], x.mut$tbl[, c(1,12, 2,3,4,5,7,8,13)])
        tbl <- tbl[order(tbl$motifName, tbl$motifRelativeScore, decreasing=TRUE),]
        tbl$signature <- sprintf("%s;%s;%s", tbl$motifName, tbl$motifStart, tbl$strand)
        tbl <- tbl [, c(1,2,10,3:9)]

        signatures.in.both <- intersect(subset(tbl, status=="mut")$signature, subset(tbl, status=="wt")$signature)
        signatures.only.in.wt <- setdiff(subset(tbl, status=="wt")$signature, subset(tbl, status=="mut")$signature)
        signatures.only.in.mut <- setdiff(subset(tbl, status=="mut")$signature, subset(tbl, status=="wt")$signature)

        tbl$assessed <- rep("failed", nrow(tbl))

        if(length(signatures.in.both) > 0) {
           indices <- sort(unlist(lapply(signatures.in.both, function(sig) grep(sig, tbl$signature))))
           tbl$assessed[indices] <- "in.both"
           }

        if(length(signatures.only.in.wt) > 0) {
           indices <- sort(unlist(lapply(signatures.only.in.wt, function(sig) grep(sig, tbl$signature))))
           tbl$assessed[indices] <- "wt.only"
           }

        if(length(signatures.only.in.mut) > 0) {
           indices <- sort(unlist(lapply(signatures.only.in.mut, function(sig) grep(sig, tbl$signature))))
           tbl$assessed[indices] <- "mut.only"
           }

        tbl$delta <- 0

           # find the mut scores for each of the "wt.only" entries, subtract from the wt score
        tbl.wt.only  <- subset(tbl, assessed=="wt.only", select=c(signature, motifRelativeScore))
        if(nrow(tbl.wt.only) > 0){
           sigs <- tbl.wt.only$signature
           tbl.mut.scores <- subset(tbl.mut.50, signature %in% sigs, select=c(signature, motifRelativeScore))
           deltas <- unlist(lapply(sigs, function(sig){wt.score  <- subset(tbl.wt.only, signature==sig)$motifRelativeScore;
                                                mut.score <- subset(tbl.mut.scores, signature==sig)$motifRelativeScore;
                                                delta <- wt.score - mut.score
                                                }))
           tbl$delta[match(sigs, tbl$signature)] <- deltas
           } # if some wt.only entries

           # find the wt scores for each of the "mut.only" entries, subtract from the mut score
        tbl.mut.only  <- subset(tbl, assessed=="mut.only", select=c(signature, motifRelativeScore))
        sigs <- tbl.mut.only$signature
        tbl.wt.scores <- subset(tbl.wt.50, signature %in% sigs, select=c(signature, motifRelativeScore))
        deltas <- unlist(lapply(sigs, function(sig){mut.score  <- subset(tbl.mut.only, signature==sig)$motifRelativeScore;
                                                    wt.score <- subset(tbl.wt.scores, signature==sig)$motifRelativeScore;
                                                    delta <- wt.score - mut.score
                                                 }))
        tbl$delta[match(sigs, tbl$signature)] <- deltas
        coi <-  c("motifName", "status", "assessed", "motifRelativeScore", "delta",
                  "signature", "chrom", "motifStart", "motifEnd", "strand",
                  "match", "tf")

        tbl <- tbl[, coi]
        tbl$variant <- variant
        tbl
        }) # assessSnp

#------------------------------------------------------------------------------------------------------------------------
