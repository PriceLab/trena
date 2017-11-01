library(trena)
print(load("../extdata/bugDemoData/solver.bug.2017-11-01.09:49:40.RData"))

solver <- EnsembleSolver(mtx, targetGene, candidates, solver.names)
run(solver)

# Error in `[.data.frame`(out.list$out.sqrtlasso, , c("beta", "gene")) :
#   undefined columns selected

# paul-shannon's sessionInfo
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Sierra 10.12.3
#
# Matrix products: default
# BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
#
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#
# attached base packages:
# [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#  [1] RUnit_0.4.31         biomaRt_2.32.1       org.Mm.eg.db_3.4.1   annotate_1.54.0      XML_3.98-1.9         AnnotationDbi_1.38.2 Biobase_2.36.2       colorspace_1.3-2     trenaViz_0.99.24     GenomicRanges_1.28.6 GenomeInfoDb_1.12.3  graph_1.54.0         Rcpp_0.12.13         BrowserViz_1.9.15    httpuv_1.3.5         jsonlite_1.5         trena_0.99.192       MotifDb_1.19.18      Biostrings_2.44.2    XVector_0.16.0       IRanges_2.10.5       S4Vectors_0.14.7     BiocGenerics_0.22.1  glmnet_2.0-13        foreach_1.4.3        Matrix_1.2-11
#
# loaded via a namespace (and not attached):
#  [1] BSgenome.Hsapiens.UCSC.hg38_1.4.1        lattice_0.20-35                          Rsamtools_1.28.0                         digest_0.6.12                            BSgenome.Mmusculus.UCSC.mm10_1.4.0       lassopv_0.1.3                            lars_1.2                                 RSQLite_2.0                              zlibbioc_1.22.0                          rlang_0.1.2                              RPostgreSQL_0.6-2                        data.table_1.10.4-2                      blob_1.1.0                               RMySQL_0.10.13                           BiocParallel_1.10.1                      SNPlocs.Hsapiens.dbSNP150.GRCh38_0.99.20 stringr_1.2.0                            igraph_1.1.2                             RCurl_1.95-4.8                           bit_1.1-12                               DelayedArray_0.2.7                       compiler_3.4.1                           rtracklayer_1.36.6                       pkgconfig_2.0.1                          vbsr_0.0.5                               SummarizedExperiment_1.6.5               tibble_1.3.4                             GenomeInfoDbData_0.99.0                  codetools_0.2-15                         matrixStats_0.52.2                       randomForest_4.6-12                      BSgenome.Hsapiens.UCSC.hg19_1.4.0        GenomicAlignments_1.12.2                 MASS_7.3-47                              bitops_1.0-6                             grid_3.4.1                               xtable_1.8-2                             DBI_0.7                                  magrittr_1.5                             stringi_1.1.5                            splitstackshape_1.4.2                    org.Hs.eg.db_3.4.1                       iterators_1.0.8                          tools_3.4.1                              bit64_0.9-7                              BSgenome_1.44.2                          flare_1.5.0                              memoise_1.1.0
