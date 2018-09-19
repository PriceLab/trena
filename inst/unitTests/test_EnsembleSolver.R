library(trena)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_EnsembleSolverConstructor()
   test_ampAD.mef2c.154tfs.278samples.ensemble()
   test_selectedSolversOnly()
   test_getSolverNames()
   test_oneSolver()
   test_invalidSolvers()

} # runTests
#----------------------------------------------------------------------------------------------------
test_EnsembleSolverConstructor <- function()
{
    printf("--- test_EnsembleSolverConstructor")

    # Construct the EnsembleSolver and check that it's correct

    #solver <- EnsembleSolver()
    mtx <- matrix(1:9,nrow=3)
    rownames(mtx) <- c("gene1","gene2","gene3")
    solver <- EnsembleSolver(mtx,targetGene = "gene1",
                             candidateRegulators = c("gene2","gene3"))
    checkEquals(class(solver)[1], "EnsembleSolver")
    checkTrue(all(c("EnsembleSolver", "Solver") %in% is(solver)))
}

# test_EnsembleSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.ensemble <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.ensemble")

   set.seed(122113)
   # Load matrix and transform via arcsinh
   load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
    mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   solver <- EnsembleSolver(mtx.asinh,target.gene,tfs)
   tbl <- run(solver)

   # Check for empirical values
   checkTrue(c("HLF") %in% tbl$gene)

} # test_ampAD.mef2c.154tfs.278samples.ensemble
#----------------------------------------------------------------------------------------------------
test_selectedSolversOnly <- function()
{
   printf("--- test_selectedSolversOnly")

   set.seed(122113)
   # Load matrix and transform via arcsinh
   load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   solvers <- c("lasso", "ridge", "lassopv", "pearson", "spearman") # "sqrtlasso",

   solver <- EnsembleSolver(mtx.asinh,target.gene,tfs,solverNames=solvers)
   tbl <- run(solver)

   # Check for empirical values
   checkTrue(c("HLF") %in% tbl$gene)

} # test_selectedSolversOnly
#----------------------------------------------------------------------------------------------------
doNot_test_pcaError <- function()
{
    printf("--- test_pcaError")

    # Take a small subset of the matrix; only 2 columns
    set.seed(122113)
    # Load matrix and transform via arcsinh
    load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))

    # Find the MEF2C row
    mef2c.idx <- which(rownames(mtx.sub) == "MEF2C")
    start.idx <- mef2c.idx - 1
    end.idx <- mef2c.idx + 1

    # Subset the matrix so it's 3 x 2
    mtx.sub <- mtx.sub[start.idx:end.idx,1:2]

    mtx.asinh <- asinh(mtx.sub)
    #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

    target.gene <- "MEF2C"
    tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
    solvers <- c("lasso", "ridge", "lassopv", "pearson", "spearman") # "sqrtlasso",

    solver <- EnsembleSolver(mtx.asinh,target.gene,tfs,solverNames=solvers)

    # Change warnings to errors
    options(warn = 2)

    checkException(run(solver), silent=TRUE)

    # Change warnings back to warnings
    options(warn = 1)

    # todo:  track down this error, add tests
    # tbl <- suppressWarnings(run(solver))
    # (03 mar 2018) fails at
    # Error in colMeans(yres) : 'x' must be an array of at least two dimensions
    #  colMeans(yres)
    #  t(t(yres) - colMeans(yres))

    # Check that pcaMax and concordance were added
    #checkTrue(ncol(tbl) == 8)
    #checkTrue(all(c("pcaMax","concordance") %in% names(tbl)))

    # Check that they're all NA
    #checkTrue(all(is.na(tbl$concordance)))
    #checkTrue(all(is.na(tbl$pcaMax)))

} # test_pcaError
#----------------------------------------------------------------------------------------------------
test_getSolverNames <- function(){

    printf("--- test_getSolverNames")

    load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
    targetGene <- "MEF2C"
    candidateRegulators <- setdiff(rownames(mtx.sub), targetGene)
    solver <- EnsembleSolver(mtx.sub, targetGene, candidateRegulators,
                             solverNames = c("lasso","randomForest"))

    solver.names <- getSolverNames(solver)

    # Test that it's what we want
    checkEquals(solver.names, c("lasso","randomForest"))
} # test_getSolverNames
#----------------------------------------------------------------------------------------------------
test_oneSolver <- function(){

    printf("--- test_oneSolver")

    load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
    targetGene <- "MEF2C"
    candidateRegulators <- setdiff(rownames(mtx.sub), targetGene)

    # Supply only a Pearson solver
    solver <- EnsembleSolver(mtx.sub, targetGene, candidateRegulators,
                             solverNames = c("pearson"), geneCutoff = 1)
    # Check for a warning
    options(warn = 2)
    checkException(run(solver), silent = TRUE)

    # Set warnings back to non-errors
    options(warn = 1)

    # Check that the output matches the Pearson output
    tbl.ens <- suppressWarnings(run(solver))

    p.solver <- PearsonSolver(mtx.sub, targetGene, candidateRegulators)
    tbl.p <- run(p.solver)

    checkEquals(names(tbl.p), names(tbl.ens))
    checkEquals(tbl.p$coefficient, tbl.ens$coefficient)
    checkEquals(rownames(tbl.p), rownames(tbl.ens))

} # test_oneSolver
#----------------------------------------------------------------------------------------------------
test_invalidSolvers <- function(){

    printf("--- test_invalidSolvers")

    load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
    targetGene <- "MEF2C"
    candidateRegulators <- setdiff(rownames(mtx.sub), targetGene)

    # Test with only an invalid solver
    solver <- EnsembleSolver(mtx.sub, targetGene, candidateRegulators,
                             solverNames = c("rudge"))

    checkException(run(solver), silent = TRUE)

    # Test with valid and invalid solvers
    options(warn = 2)
    solver <- EnsembleSolver(mtx.sub, targetGene, candidateRegulators,
                             solverNames = c("lasso","ridge","parson"))
    checkException(run(solver), silent = TRUE)

    # Test to make sure the output is just lasso and ridge
    options(warn = 1)
    set.seed(11451)
    tbl <- suppressWarnings(run(solver))
    checkEquals(colnames(tbl), c("gene", "betaLasso", "betaRidge"))
    checkTrue(!("parson" %in% colnames(tbl)))

} # test_invalidSolvers
#----------------------------------------------------------------------------------------------------
# pshannon (19 sep 2018): scoring multiple solvers together is no longer done
#                         within the EnsemblSolver class: it is a post-processing step,
#                         of some complexity, with no single analytical "right way'
#                         see the soon-to-appear class to be named something like "EnsembleScorer"
# pshannon (16 jun 2018):
#   when cory uses the new trenaSGM package, his model shows strong disagreement between
#   rfScore and pcaMax.  our experience has been that these two scores roughly agree.
#   here we explore, determine and fix that error, using some simple precalculated
#   ensembleSolver gene model results
#
#   3 serialize datasets are used here:
#      1) pcaMaxBug.01.RData: tbl.all & tbl.scale, illustrates the rfScore/pcaMax disagreement
#      2) pcaMaxBug.02.RData: tbl.all & tbl.scale, rfScore & pcaMax roughly agree
#      3) pacaMaxBut.03.Rdata: much larger model because TFClass mapping is used,
#         here too rfScore & pcaMax
#
doNot_test_.addEnsembleScore <- function()
{
       # demonstrate the problem
   load(system.file(package="trena", "extdata", "pcaMaxBug.01.RData"))
   tbl.new <- trena:::.addEnsembleScores(tbl.scale, tbl.all)
   genes.by.rf <- tbl.new[order(tbl.new$rfScore, decreasing=TRUE), "gene"][1:6]
   genes.by.pcaMax <- tbl.new[order(tbl.new$pcaMax, decreasing=TRUE), "gene"][1:6]
   relative.ordering <- match(genes.by.rf, genes.by.pcaMax)  # 2 4 5 6 1 3

      # is scaling the problem?   try it here
   mtx.scale <- as.matrix(tbl.all[, -1])   # remove the gene column
   rownames(mtx.scale) <- tbl.all$gene

   scale.scores.normal.distribution <- function(vec){
      vec.center <- median(vec)
      vec.scale  <- mad(vec)
      scale(vec, center=vec.center, scale=vec.scale)
      }

   scale.scores.poisson.distribution <- function(vec){
      vec.center <- median(vec)
      vec.scale <- sqrt(mean(vec*vec))
      scale(vec, center=vec.center, scale=vec.scale)
      }

   mtx.scale[, "betaLasso"] <- scale.scores.normal.distribution(mtx.scale[, "betaLasso"])
   mtx.scale[, "rfScore"] <- scale.scores.poisson.distribution(mtx.scale[, "rfScore"])
   mtx.scale[, "betaLasso"] <- scale.scores.normal.distribution(mtx.scale[,"betaLasso"])
   mtx.scale[, "lassoPValue"] <- scale.scores.normal.distribution(-log10(mtx.scale[,"lassoPValue"]))
   mtx.scale[, "pearsonCoeff"] <- scale.scores.normal.distribution(mtx.scale[,"pearsonCoeff"])
   mtx.scale[, "betaSqrtLasso"] <- scale.scores.normal.distribution(mtx.scale[,"betaSqrtLasso"])
   mtx.scale[, "rfScore"] <- scale.scores.normal.distribution(mtx.scale[,"rfScore"])
   mtx.scale[, "betaRidge"] <- scale.scores.normal.distribution(mtx.scale[,"betaRidge"])
   mtx.scale[, "spearmanCoeff"] <- scale.scores.normal.distribution(mtx.scale[,"spearmanCoeff"])

    pca <- stats::prcomp(mtx.scale, center=FALSE, scale.=FALSE)

    pca$x <- pca$x / sqrt(length(which(pca$sdev > 0.1)))
    pcaMax <- apply(pca$x[, pca$sdev > 0.1, drop=FALSE],1, function(x) {sqrt(mean(x*x))})
    pcaMax <- as.data.frame(pcaMax)
    pcaMax$gene <- rownames(pcaMax)

       #         pcaMax  gene
       # CEBPA 0.6601057 CEBPA
       # IKZF1 2.3129081 IKZF1
       # IRF2  0.5121313  IRF2
       # IRF8  0.2370583  IRF8
       # NR6A1 2.9354326 NR6A1
       # TAL1  0.1831896  TAL1


} # test_.addEnsembleScore
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
