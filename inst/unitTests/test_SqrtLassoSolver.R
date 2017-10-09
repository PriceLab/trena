library(trena)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_SqrtLassoSolverConstructor()
   test_ampAD.mef2c.154tfs.278samples.sqrtlasso()
   test_lambda.sqrtlasso()
   test_nCores.sqrtlasso()
   
} # runTests
#----------------------------------------------------------------------------------------------------
test_SqrtLassoSolverConstructor <- function()
{
    printf("--- test_SqrtLassoSolverConstructor")

    mtx <- matrix(rnorm(100),nrow=10)   
    rownames(mtx) <- paste0("gene",1:10)
    solver <- SqrtLassoSolver(mtx,targetGene = "gene1",                          
                          candidateRegulators = c("gene2","gene3"))    
   
    checkEquals(class(solver)[1], "SqrtLassoSolver")    
    checkTrue(all(c("SqrtLassoSolver", "Solver") %in% is(solver)))
}

# test_SqrtLassoSolverConstructor   
#----------------------------------------------------------------------------------------------------    
test_ampAD.mef2c.154tfs.278samples.sqrtlasso <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.sqrtlasso")

   # Load matrix and transform via arcsinh
   load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))   
   mtx.asinh <- asinh(mtx.sub)
   target.gene <- "MEF2C"
   # Use only 30 genes
   tfs <- setdiff(rownames(mtx.asinh)[1:10], "MEF2C")
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   set.seed(10)
   sqrt.solver <- SqrtLassoSolver(mtx.asinh, target.gene, tfs)
   tbl <- run(sqrt.solver)

   # Check for empirical values
   tbl <- tbl[order(abs(tbl$beta), decreasing=TRUE),, drop = FALSE]
   expected.genes <- sort(c("ATF2","ADNP2","ARX","ATF7"))
   actual.genes <- sort(rownames(subset(tbl, abs(beta) > 0.2)))   
   checkTrue(all(expected.genes %in% actual.genes))
   checkTrue(max(tbl$beta) < 1.45)
   checkTrue(min(tbl$beta) > -0.25)
   
} # test_ampAD.mef2c.154tfs.278samples.sqrtlasso
#----------------------------------------------------------------------------------------------------
test_lambda.sqrtlasso <- function()
{
   printf("--- test_lambda.sqrtlasso")

   # Load matrix and transform via arcsinh
   load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))   
   mtx.asinh <- asinh(mtx.sub)
   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh)[1:10], "MEF2C")
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   set.seed(1144)
   sqrt.solver <- SqrtLassoSolver(mtx.asinh, target.gene, tfs, lambda = 0.1)
   tbl <- run(sqrt.solver)   
   
   # Check for empirical values
   tbl <- tbl[order(abs(tbl$beta), decreasing=TRUE),, drop = FALSE]
   expected.genes <- sort(c("ADNP","ADNP2","ARX","ATF2", "ATF7", "BCL6"))
   actual.genes <- sort(rownames(subset(tbl, abs(beta) > 0.1)))
   checkEquals(expected.genes,actual.genes)
   checkTrue(max(tbl$beta) < 1.45)
   checkTrue(min(tbl$beta) > -0.25)
   
} # test_lambda.sqrtlasso
#----------------------------------------------------------------------------------------------------
test_nCores.sqrtlasso <- function()
{
   printf("--- test_nCores.sqrtlasso")

   # Load matrix and transform via arcsinh
   load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))   
   mtx.asinh <- asinh(mtx.sub)
   target.gene <- "MEF2C"
   # Use only 30 genes
   tfs <- setdiff(rownames(mtx.asinh)[1:30], "MEF2C")
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   set.seed(10)
   sqrt.solver <- SqrtLassoSolver(mtx.asinh, target.gene, tfs, nCores = parallel::detectCores() - 1)
   tbl <- run(sqrt.solver)

   # Check for empirical values
   tbl <- tbl[order(abs(tbl$beta), decreasing=TRUE),, drop = FALSE]
   expected.genes <- sort(c("ATF2","CUX1","ESRRG","FOXD4L1"))
   actual.genes <- sort(rownames(subset(tbl, abs(beta) > 0.2)))   
   checkTrue(all(expected.genes %in% actual.genes))
   checkTrue(max(tbl$beta) < 0.75)
   checkTrue(min(tbl$beta) > -0.55)
   
} # test_ncores.sqrtlasso
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
