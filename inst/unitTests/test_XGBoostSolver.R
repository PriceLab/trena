library(trena)
library(RUnit)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_XGBoostSolverConstructor()
   test_ampAD.mef2c.154tfs.278samples.XGBoost()

} # runTests
#----------------------------------------------------------------------------------------------------
test_XGBoostSolverConstructor <- function()
{
    message(sprintf("--- test_XGBoostSolverConstructor"))

    mtx <- matrix(1:9,nrow=3)
    rownames(mtx) <- c("gene1","gene2","gene3")
    solver <- XGBoostSolver(mtx,targetGene = "gene1",
                            candidateRegulators = c("gene2","gene3"))

    checkEquals(class(solver)[1], "XGBoostSolver")
    checkTrue(all(c("XGBoostSolver", "Solver") %in% is(solver)))

} # test_XGBoostSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.XGBoost <- function()
{
   message(sprintf("--- test_ampAD.mef2c.154tfs.278samples.XGBoost"))

   # Load matrix and transform via arcsinh
   load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))

   mtx.asinh <- asinh(mtx.sub)
   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   XGBoost.solver <- XGBoostSolver(mtx.asinh, target.gene, tfs)
   tbl <- run(XGBoost.solver)

   # Check for empirical values
   checkTrue(nrow(subset(tbl, abs(Importance) > 0.013)) == 9)

} # test_ampAD.mef2c.154tfs.278samples.XGBoost
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
