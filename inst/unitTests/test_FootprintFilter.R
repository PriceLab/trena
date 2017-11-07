library(trena)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_FootprintFilter.byRegion()
    test_FootprintFilter.byTwoRegions()
    
} # runTests
#----------------------------------------------------------------------------------------------------
# should exactly replicate the geneCentered results above
test_FootprintFilter.byRegion <- function()
{
    printf("--- test_FootprintFilter.byRegion")
    
    db.address <- system.file(package="trena", "extdata")
    genome.db.uri    <- paste("sqlite:/", db.address, "mef2c.neighborhood.hg38.gtfAnnotation.db",  sep = "/")
    footprint.db.uri <- paste("sqlite:/", db.address, "mef2c.neigborhood.hg38.footprints.db", sep = "/")
    #target.gene <- "MEF2C"
    
    mef2c.tss <- 88904257   # minus strand
    chrom <- "chr5"
    start <- mef2c.tss - 1000
    end   <- mef2c.tss + 1000
    tbl.regions <- data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE)
    
    recipe <- list(genomeDB=genome.db.uri,
                   footprintDB=footprint.db.uri,
                   regions=tbl.regions)
    
    filter <- FootprintFilter(genome.db.uri,
                              footprint.db.uri,
                              tbl.regions,
                              quiet=TRUE)
    
    tbl.out <- getCandidates(filter)[[1]]
    min.fp.start <- min(tbl.out$start)
    max.fp.end   <- max(tbl.out$end)
    checkTrue(min.fp.start > (mef2c.tss - 1000))
    checkTrue(max.fp.end <  (mef2c.tss + 1000))
    span <- max.fp.end - min.fp.start   # less than 2000, more than nothing
    checkTrue(span > 100)
    checkTrue(span <= 20001)
    
} # test_FootprintFilter.byRegion
#----------------------------------------------------------------------------------------------------
# get footprints for tss +/- 1000, by specifying two contiguous regions, one 500 bp, one 1500 bp
test_FootprintFilter.byTwoRegions <- function()
{
    printf("--- test_FootprintFilter.byTwoRegions")
    
    db.address <- system.file(package="trena", "extdata")
    genome.db.uri    <- paste("sqlite:/", db.address, "mef2c.neighborhood.hg38.gtfAnnotation.db",  sep = "/")
    footprint.db.uri <- paste("sqlite:/", db.address, "mef2c.neigborhood.hg38.footprints.db", sep = "/")
    
    mef2c.tss <- 88904257   # minus strand
    
    tbl.regions <- data.frame(chrom=c("chr5", "chr5"),
                              start=c(mef2c.tss - 1000, mef2c.tss - 500),
                              end=c(mef2c.tss - 500, mef2c.tss + 1000),
                              stringsAsFactors=FALSE)
    
    recipe <- list(genomeDB=genome.db.uri,
                   footprintDB=footprint.db.uri,
                   regions=tbl.regions)
    
    filter <- FootprintFilter(recipe$genomeDB,
                              recipe$footprintDB,
                              recipe$regions,
                              quiet=TRUE)
    
    list.out <- getCandidates(filter)
    checkTrue(length(list.out) == 2)
    tbl.out <- rbind(list.out[[1]],list.out[[2]])
    min.fp.start <- min(tbl.out$fp_start)
    max.fp.end   <- max(tbl.out$fp_end)
    checkTrue(min.fp.start > (mef2c.tss - 1020))
    checkTrue(max.fp.end <  (mef2c.tss + 1000))
    span <- max.fp.end - min.fp.start   # less than 2000, more than nothing
    checkTrue(span > 100)
    checkTrue(span <= 20001)
    
    checkTrue(nrow(tbl.out) > 130)
    checkTrue(nrow(tbl.out) < 150)
    
} # test_FootprintFilter.byTwoRegions
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
