printf <- function(...) print(noquote(sprintf(...)))

install.packages("BiocManager", repos="https://cran.rstudio.com")
biocGet <- function(pkgs){
   library(BiocManager)
   BiocManager::install(pkgs, update=TRUE, ask=FALSE, type="source", version = "3.16", force=TRUE)
   }


code.pkgs <- c(
    "AnnotationDbi",
    "BiocGenerics",
    "BiocParallel",
    "biomaRt",
    "Biostrings",
    "BSgenome",
    "BSgenome.Hsapiens.UCSC.hg19",
    "BSgenome.Hsapiens.UCSC.hg38",
    "BSgenome.Mmusculus.UCSC.mm10",
    "data.table",
    "DBI",
    "formatR",
    "GenomicRanges",
    "glmnet",
    "knitr",
    "lassopv",
    "markdown",
    "methods",
    "MotifDb",
    "org.Hs.eg.db",
    "plyr",
    "R6",
    "randomForest",
    "rmarkdown",
    "RMySQL",
    "RPostgreSQL",
    "RSQLite",
    "RUnit",
    "SNPlocs.Hsapiens.dbSNP150.GRCh38",
    "utils",
    "WGCNA",
    "xgboost")



for(code.pkg in code.pkgs){
   suppressWarnings(
      needed <- !require(code.pkg, character.only=TRUE, lib.loc=my.user.library, quiet=TRUE)
      )
   printf("%s needed? %s", code.pkg, needed)
   if(needed)
      biocGet(code.pkg)
   } # for



