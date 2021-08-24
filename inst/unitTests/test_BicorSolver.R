library(trena)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_BicorSolverConstructor()
   test_ampAD.mef2c.154tfs.278samples.bicor()

} # runTests
#----------------------------------------------------------------------------------------------------
test_BicorSolverConstructor <- function()
{
    printf("--- test_BicorSolverConstructor")

    mtx <- matrix(1:9,nrow=3)
    rownames(mtx) <- c("gene1","gene2","gene3")
    solver <- BicorSolver(mtx,targetGene = "gene1", candidateRegulators = c("gene2","gene3"))

    checkEquals(class(solver)[1], "BicorSolver")
    checkTrue(all(c("BicorSolver", "Solver") %in% is(solver)))

} # test_BicorSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.bicor <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.bicor")

   # Load matrix and transform via arcsinh
   load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))

   mtx.asinh <- asinh(mtx.sub)
   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   bicor.solver <- BicorSolver(mtx.asinh, target.gene, tfs)
   tbl <- run(bicor.solver)

   # Check for empirical values
   checkTrue(nrow(subset(tbl, abs(coefficient) > 0.8)) == 9)

} # test_ampAD.mef2c.154tfs.278samples.bicor
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
