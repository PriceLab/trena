library(trena)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_CandidateFilter()
   
} # runTests
#----------------------------------------------------------------------------------------------------
test_CandidateFilter <- function()
{
    printf("--- test_CandidateFilter")

    # Create an empty filter
    candidate.filter <- CandidateFilter(quiet = TRUE)
    checkEquals(is(candidate.filter), "CandidateFilter")

} #test_CandidateFilter
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
