biocGet <- function(pkgs){
   library(BiocManager)
   BiocManager::install(pkgs)
   }

printf <- function(...) print(noquote(sprintf(...)))
code.pkgs <- c("GenomicRanges",
               "Biostrings",
               "jsonlite",
               "httpuv",
               "BiocGenerics",
               "BiocParallel",
               "DBI",
               "RPostgreSQL",
               "RMySQL",
               "RSQLite",
               "glmnet",
               "lassopv",
               "randomForest",
               "flare",
               "vbsr",
               "stringr",
               "httpuv",
               "colorspace",
               "annotate",
               "shiny",
               "DT",
               "htmlwidgets",
               "later",
               "splitstackshape",
               "RUnit")

for(code.pkg in code.pkgs){
   suppressWarnings(
      needed <- !require(code.pkg, character.only=TRUE)
      )
   printf("%s needed? %s", code.pkg, needed)
   if(needed)
      biocGet(code.pkg)
   } # for
