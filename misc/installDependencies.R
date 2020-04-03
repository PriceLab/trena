biocGet <- function(pkgs){
   library(BiocManager)
   BiocManager::install(pkgs)
   }

printf <- function(...) print(noquote(sprintf(...)))
code.pkgs <- c(
    "annotate",
    "AnnotationDb",
    "BiocGenerics",
    "BiocParallel",
    "biomaRt",
    "Biostrings",
    "BSgenome",
    "BSgenome.Athaliana.TAIR.TAIR9",
    "BSgenome.Hsapiens.UCSC.hg19",
    "BSgenome.Hsapiens.UCSC.hg38",
    "BSgenome.Mmusculus.UCSC.mm10",
    "BSgenome.Scerevisiae.UCSC.sacCer3",
    "colorspace",
    "DBI",
    "DT",
    "flare",
    "GenomicRanges",
    "glmnet",
    "glmnet (>= 2.0.3)",
    "htmlwidgets",
    "httpuv",
    "jsonlite",
    "knitr",
    "lassopv",
    "later",
    "org.Hs.eg.db",
    "org.Mm.eg.db",
    "plyr",
    "randomForest",
    "rmarkdown",
    "RMySQL",
    "RPostgreSQL",
    "RSQLite",
    "RUnit",
    "shiny",
    "splitstackshape",
    "stringr",
    "vbsr",
    "xgboost")

for(code.pkg in code.pkgs){
   suppressWarnings(
      needed <- !require(code.pkg, character.only=TRUE)
      )
   printf("%s needed? %s", code.pkg, needed)
   if(needed)
      biocGet(code.pkg)
   } # for
