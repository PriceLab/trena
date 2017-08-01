library(trena)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_VarianceFilter()

} # runTests
#----------------------------------------------------------------------------------------------------
test_VarianceFilter <- function()
{
    printf("--- test_VarianceFilter")

    # Create dummy data and filter it based on variance
    load(system.file(package="trena", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
    var.filter <- VarianceFilter(mtx.assay = mtx.sub,
                                 targetGene ="MEF2C",
                                 varSize = 0.5)
    tf.list <- getCandidates(var.filter)
    
    checkTrue(length(tf.list$tfs) == 1)

} # test_VarianceFilter
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
