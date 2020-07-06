#----------------------------------------------------------------------------------------------------
#' Parse a string containing a chromosome and location on the genome
#'
#' Given a string of the format "chromosome:start-end", pull out the three values contained in
#' the string and return them as a named list
#'
#' @param chromLocString A string of the format "chromosome:start-end"
#'
#' @return A named list containing the chromosome, start, and end from the supplied string
#'
#' @export
#'
#' @examples
#' # Note: Both examples return :
#' # list(chrom="chr10", start=118441047, end=118445892)
#'
#' # Parse a string containing the "chr" prefix for chromosome
#' chrom.list <- parseChromLocString("chr10:118441047-118445892")
#'
#' # Parse a string without the "chr" prefix in the chromosome
#' chrom.list <- parseChromLocString("10:118441047-118445892")
#'
#' @aliases parseChromLocString
#' @rdname parseChromLocString

parseChromLocString <- function(chromLocString)
{
    chromLocString <- gsub(",", "", chromLocString);
    tokens.0 <- strsplit(chromLocString, ":", fixed=TRUE)[[1]]
    stopifnot(length(tokens.0) == 2)
    chrom <- tokens.0[1]
    if(!grepl("chr", chrom))
        chrom <- sprintf("chr%s", chrom)

    tokens.1 <- strsplit(tokens.0[2], "-")[[1]]
    stopifnot(length(tokens.1) == 2)
    start <- as.integer(tokens.1[1])
    end <- as.integer(tokens.1[2])

    return(list(chrom=chrom, start=start, end=end))

} # parseChromLocString
#----------------------------------------------------------------------------------------------------
#' Parse a string containing the information for connecting to a database
#'
#' Given a string of the format "database_driver://host/database_name", pull out the 3 pieces of
#' information and return them as a named list.
#'
#' @param database.uri A string of the format "database_driver://host/database_name"
#'
#' @return A named list containing the driver, host, and database name from the
#' supplied string
#'
#' @export
#'
#' @examples
#'
#' # Parse a URI for a local PostgreSQL database called "gtf"
#' parsed.URI <- parseDatabaseUri("postgres://localhost/gtf")
#'
#' # Parse a URI for the included SQLite database "vrk2.neighborhood.hg38.gtfAnnotation.db" in the local documents folder
#' db.address <- system.file(package="trena", "extdata")
#' genome.db.uri    <- paste("sqlite:/", db.address, "vrk2.neighborhood.hg38.gtfAnnotation.db",  sep = "/")
#' parsed.URI <- parseDatabaseUri(genome.db.uri)
#'
#' @aliases parseDatabaseUri
#' @rdname parseDatabaseUri

parseDatabaseUri <- function(database.uri)
{
    topLevel.tokens <- strsplit(database.uri, "://")[[1]]
    database.brand <- topLevel.tokens[1]
    secondLevel.tokens <- strsplit(topLevel.tokens[2], "/(?=[^/]+$)", perl = TRUE)[[1]]
    host <- secondLevel.tokens[1]
    database.name <- secondLevel.tokens[2]

    list(brand=database.brand, host=host, name=database.name)

} # parseDatabaseUri
#----------------------------------------------------------------------------------------------------
#' Get the available solvers for use in trena
#'
#' Retrieve the vector of different methods that can be used as solvers in trena. Solver names in
#' the returned vector correspond to the exact capitalization used in their Solver subclasses (i.e.
#' a Solver object using LASSO is LassoSolver)
#'
#' @return A character vector of all solvers currently available in trena.
#'
#' @export
#'
#' @examples
#' all.solvers <- getAvailableSolvers()
#'
#' @aliases getAvailableSolvers
#' @rdname getAvailableSolvers

getAvailableSolvers <- function() {

    availableSolvers <- c("BayesSpike", "LassoPV", "Lasso", "Pearson",
                          "RandomForest", "Ridge", "Spearman")

    return(availableSolvers)

} # getAvailableSolvers
#----------------------------------------------------------------------------------------------------
