library(trena)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
if(!exists("tbl.mef2c")){
   print(load(system.file(package="trena", "extdata", "tbl.mef2c.model.RData")))
   }
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_.normalize_randomForest()
   test_.normalize_pvals()
   test_.normalize_betas()
   test_.normalize_correlations()
   test_normalizeModel()
   test_addStats()

} # runTests
#----------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")
   pcaMaxCalculator <- PCAMax(tbl.mef2c)
   checkEquals(is(pcaMaxCalculator), "PCAMax")

} # test_constructor
#----------------------------------------------------------------------------------------------------
test_.normalize_randomForest <- function()
{
   printf("--- test_normalize_randomFores")
   desired.max <- 10   # for random forest, the 3rd quartile is scaled up or down to 75% of this
   vn <- trena:::.normalize_randomForest(tbl.mef2c$rfScore, normalizing.max=desired.max)
   checkEqualsNumeric(fivenum(vn)[4], 0.75 * desired.max)

} # test_.mormalize_randomForest
#----------------------------------------------------------------------------------------------------
test_.normalize_pvals <- function()
{
   printf("--- test_normalize_pvals")
   desired.max <- 10
   vn <- trena:::.normalize_pval(tbl.mef2c$lassoPValue, normalizing.max=desired.max)
   checkEqualsNumeric(max(vn), desired.max)

} # test_.mormalize_pvals
#----------------------------------------------------------------------------------------------------
test_.normalize_betas <- function()
{
   printf("--- test_normalize_betas")
   vn.0 <- trena:::.normalize_betaValues(tbl.mef2c$betaLasso, normalizing.max=10)
   vn.1 <- trena:::.normalize_betaValues(tbl.mef2c$betaRidge, normalizing.max=10)
   vn.2 <- trena:::.normalize_betaValues(tbl.mef2c$betaSqrtLasso, normalizing.max=10)

   checkEquals(as.numeric(lapply(list(vn.0, vn.1, vn.2), max)), c(10, 10, 10))

} # test_.normalize_betas
#----------------------------------------------------------------------------------------------------
test_.normalize_correlations <- function()
{
   printf("--- test_normalize_betas")
   vn.0 <- trena:::.normalize_correlationValues(tbl.mef2c$spearmanCoeff, normalizing.max=10)
   vn.1 <- trena:::.normalize_correlationValues(tbl.mef2c$pearsonCoeff, normalizing.max=10)
   checkTrue(max(vn.0) > 8)
   checkTrue(max(vn.1) > 8)

} # test_.normalize_correlations
#----------------------------------------------------------------------------------------------------
test_pca <- function()
{


} # test_pca
#----------------------------------------------------------------------------------------------------
test_normalizeModel <- function()
{
   printf("--- test_normalizeModel")

   x <- PCAMax(tbl.mef2c)
   mtx <- normalizeModel(x, normalizing.max=10)
   checkEquals(dim(mtx), c(20,7))
   coi <- c("betaLasso", "lassoPValue", "pearsonCoeff", "rfScore",
            "betaRidge", "spearmanCoeff", "betaSqrtLasso")
   checkEquals(colnames(mtx), coi)
   mtx.summary <- apply(mtx, 2, fivenum)

} # test_normalizeModel
#----------------------------------------------------------------------------------------------------
test_addStats <- function()
{
   printf("--- test_addStats")

   x <- PCAMax(tbl.mef2c)
   mtx <- normalizeModel(x, normalizing.max=10)
   tbl.05 <- addStats(x, varianceToInclude=0.5, scalePCA=TRUE)
   tbl.10 <- addStats(x, varianceToInclude=.99, scalePCA=FALSE)
   checkEquals(head(tbl.05$tf.hgnc), c("GABPA", "SMAD5", "STAT4", "TCF12", "TBR1", "PKNOX2"))
   checkEquals(head(tbl.10$tf.hgnc), c("GABPA", "SMAD5", "TCF12", "STAT4", "TBR1", "PKNOX2"))

   # mtx.summary <- apply(mtx, 2, fivenum)

} # test_normalizeModel
#----------------------------------------------------------------------------------------------------

