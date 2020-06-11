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
            "betaRidge", "spearmanCoeff")
   checkEquals(colnames(mtx), coi)
   mtx.summary <- apply(mtx, 2, fivenum)

} # test_normalizeModel
#----------------------------------------------------------------------------------------------------
test_addStats <- function()
{
   printf("--- test_addStats")

   coi <- c(2,5:11,15:17)
   coi.2 <- c(2,5:11,15)
   x <- PCAMax(tbl.mef2c)
   mtx <- normalizeModel(x, normalizing.max=10)
   tbl.05 <- addStatsSimple(x, varianceToInclude=0.5, scalePCA=TRUE, quiet=FALSE)
   tbl.05 <- tbl.05[, coi]
   tbl.05 <- tbl.05[order(tbl.05$score, decreasing=TRUE),]
   tbl.05.2 <- addStats(x, varianceToInclude=0.5, scalePCA=TRUE)
   tbl.10 <- addStats(x, varianceToInclude=.99, scalePCA=FALSE)
   checkEquals(head(tbl.05$tf.hgnc), c("GABPA", "SMAD5", "STAT4", "TCF12", "TBR1", "PKNOX2"))
   checkEquals(head(tbl.10$tf.hgnc), c("GABPA", "SMAD5", "TCF12", "STAT4", "TBR1", "PKNOX2"))

   # mtx.summary <- apply(mtx, 2, fivenum)

} # test_normalizeModel
#----------------------------------------------------------------------------------------------------
explore <- function ()
{
   print(load("~/github/eqtlTrenaNotebooks/mef2c/trena/data/mtx.withDimers.cer.ros.tcx.RData"))
   tfs <- subset(tbl.mef2c, abs(pearsonCoeff) > 0.7)$tf.hgnc
   noquote(intersect(tfs, rownames(mtx.tcx)))
   fit <- lm(formula = MEF2C ~ STAT4 + PKNOX2 + HLF + GABPA + SMAD5 + TCF12 + SATB2 + EMX1 + 0,
   #fit <- lm(formula = MEF2C ~ STAT4 + TCF12 + GABPA + PKNOX2,
    fit <- lm(formula = MEF2C~STAT4+TBR1+PKNOX2+HLF+EGR3+GABPA+PAX7+SMAD5+TCF12+SATB2+SATB1+EMX1+MKX+TSHZ3+MEF2D+DLX5,
          data=as.data.frame(t(mtx.tcx)))
   tbl.coef <- as.data.frame(summary(fit)$coefficients)
   colnames(tbl.coef) <- c("estimate", "error", "t", "pval")
   tbl.coef <- tbl.coef[order(tbl.coef$pval, decreasing=FALSE),]

} # explore
#----------------------------------------------------------------------------------------------------
simulatedData <- function()
{
   printf("--- simulatedData")
   set.seed(37)
   mtx <- matrix(nrow=3, ncol=100, dimnames=list(paste0("g", 1:3), paste0("S", 1:100)))
   mtx[1,] <- round(runif(100, -3, 3), 2)
   mtx[2,] <- jitter(mtx[1,], amount=0.11)
   mtx[3,] <- jitter(mtx[2,], amount=0.1)
   cor(mtx[1,], mtx[2,])
   cor(mtx[1,], mtx[3,])
   cor(mtx[2,], mtx[3,])

      #------------------------------------------------------------
      # linear regression with base R's lm
      #------------------------------------------------------------

   fit <- lm(formula = g1 ~ g2 + g3 + 0, data=as.data.frame(t(mtx)))
   tbl.coef <- as.data.frame(summary(fit)$coefficients)
   colnames(tbl.coef) <- c("estimate", "error", "t", "pval")
   tbl.coef
      #      estimate     error         t         pval
      # g2 0.95095814 0.1232547 7.7153921 1.021649e-11
      # g3 0.04690956 0.1229210 0.3816237 7.035667e-01

      # with some earlier matrix, got these results (closer to what I expected)
      #     estimate      error        t         pval
      # g2 0.5045317 0.05157506 9.782474 3.589792e-16
      # g3 0.4896006 0.05161448 9.485722 1.582422e-15

      #------------------------------------------------------------
      # use our wrapper around lasso.  ridge (get all predictors) with alpha=0
      #------------------------------------------------------------

   lasso <- LassoSolver(mtx, "g1", c("g2", "g3"), alpha=0, lambda=1)
   run(lasso)
      #         beta   intercept
      # g2 0.5084780 0.001498278
      # g3 0.4710166 0.001498278

    rf <- RandomForestSolver(mtx,targetGene = "g1", candidateRegulators = c("g2","g3"))
    run(rf)
      #    IncNodePurity
      # g3      138.4487
      # g2      131.0454

} # simulatedData
#----------------------------------------------------------------------------------------------------
mef2c.and.tcx <- function()
{
   if(!exists("mtx.tcx"))
      load("~/github/eqtlTrenaNotebooks/mef2c/trena/data/mtx.withDimers.cer.ros.tcx.RData")
   tfs <- subset(tbl.mef2c, abs(pearsonCoeff) > 0.7)$tf.hgnc
   noquote(intersect(tfs, rownames(mtx.tcx)))
   #fit <- lm(formula = MEF2C ~ STAT4 + PKNOX2 + HLF + GABPA + SMAD5 + TCF12 + SATB2 + EMX1 + 0,
   #fit <- lm(formula = MEF2C ~ STAT4 + TCF12 + GABPA + PKNOX2,
   fit <- lm(formula = MEF2C~STAT4+TBR1+PKNOX2+HLF+EGR3+GABPA+PAX7+SMAD5+TCF12+SATB2+SATB1+EMX1+MKX+TSHZ3+MEF2D+DLX5,
             data=as.data.frame(t(mtx.tcx)))
   tbl.lm.coef <- as.data.frame(summary(fit)$coefficients)
   colnames(tbl.lm.coef) <- c("estimate", "error", "t", "pval")
   tbl.lm.coef <- round(tbl.lm.coef[order(tbl.lm.coef$pval, decreasing=FALSE),],2)
   tbl.lm.coef
   lm.scaler <- 100/max(tbl.lm.coef$estimate)
   tbl.lm.coef$scaled <- round(lm.scaler * tbl.lm.coef$estimate,2)
   tbl.lm.coef <- tbl.lm.coef[order(abs(tbl.lm.coef$scaled), decreasing=TRUE),]
      # get rid of (Intercept) row
   deleter <- grep("Intercept", rownames(tbl.lm.coef))
   if(length(deleter) > 0)
      tbl.lm.coef <- tbl.lm.coef[-deleter,]

   genes <- c("MEF2C", tfs)
   setdiff(genes, rownames(mtx.tcx))
   lasso <- LassoSolver(mtx.tcx[genes,], "MEF2C", tfs, alpha=0, lambda=1)
   tbl.lasso.coef <- round(run(lasso),2)
   tbl.lasso.coef
   lasso.factor <- 100/max(abs(tbl.lasso.coef$beta))
   tbl.lasso.coef$scaled <- round(lasso.factor * tbl.lasso.coef$beta, 2)
   tbl.lasso.coef <- tbl.lasso.coef[order(abs(tbl.lasso.coef$scaled), decreasing=TRUE),]

   tbl.combined <- cbind(tbl.lm.coef, tbl.lasso.coef[rownames(tbl.lm.coef),])[, c(5, 8)]
   colnames(tbl.combined) <- c("lm", "lasso")

     #           lm  lasso
     # PKNOX2 100.00  77.78
     # TBR1   100.00  88.89
     # STAT4   83.33 100.00
     # EGR3    66.67  77.78
     # HLF     66.67  77.78
     # TCF12  -61.11 -88.89
     # SMAD5  -55.56 -88.89
     # GABPA  -44.44 -88.89
     # PAX7    38.89  66.67
     # SATB1  -33.33  33.33
     # EMX1   -33.33  44.44
     # DLX5    27.78  77.78
     # MKX     22.22  44.44
     # TSHZ3    5.56  44.44
     # SATB2    5.56  55.56
     # MEF2D    0.00  55.56


indices <- match(rownames(tbl.lasso.coef), rownames(tbl.lm.coef))

} # mef2c.and.tcx
#----------------------------------------------------------------------------------------------------

