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
   test_pcaError()
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
   checkTrue(min(tbl$pcaMax) > 0.6)
   checkTrue(max(tbl$pcaMax) < 2.5)
   checkTrue(min(tbl$concordance) > 0.3)
   checkTrue(max(tbl$concordance) < 0.55)
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
   checkTrue(min(tbl$pcaMax) > 0.8)
   checkTrue(max(tbl$pcaMax) < 1.9)
   checkTrue(min(tbl$concordance) > 0.35)
   checkTrue(max(tbl$concordance) < 0.55)
   checkTrue(c("HLF") %in% tbl$gene)

} # test_selectedSolversOnly
#----------------------------------------------------------------------------------------------------
test_pcaError <- function()
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
    
    checkException(run(solver), silent =TRUE)

    # Change warnings back to warnings
    options(warn = 1)
    
    tbl <- suppressWarnings(run(solver))
    
    # Check that pcaMax and concordance were added
    checkTrue(ncol(tbl) == 8)
    checkTrue(all(c("pcaMax","concordance") %in% names(tbl)))

    # Check that they're all NA
    checkTrue(all(is.na(tbl$concordance)))
    checkTrue(all(is.na(tbl$pcaMax)))                                     

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
    checkEquals(names(tbl), c("gene", "betaLasso", "betaRidge", "concordance", "pcaMax"))
    checkTrue(max(tbl$pcaMax) > 22.5)
    checkTrue(max(tbl$concordance) > 0.7)
    checkTrue(!("parson" %in% names(tbl)))
        
} # test_invalidSolvers
#----------------------------------------------------------------------------------------------------
    
if(!interactive()) runTests()
