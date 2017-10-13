# Load up the relevant libraries
library(trena)
library(RPostgreSQL)

# Create a function that runs on either BDDS and Khaleesi, pulling times for each

timeFpFilter <- function(host){
    
    # Set the parameters for whichever database is provided    
    genome.db.uri <- paste0("postgres://",host,"/hg38")
    footprint.db.uri <- paste0("postgres://", host, "/brain_hint_20")
    chrom <- "chr5"
    start <- 88903373
    end   <- 88907372
    
    # Create a FootprintFilter
    tbl.regions <- data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE)
    recipe <- list(genomeDB=genome.db.uri,
                   footprintDB=footprint.db.uri,
                   regions=tbl.regions)
    
    fp.filter <- FootprintFilter(genome.db.uri,
                                 footprint.db.uri,
                                 tbl.regions,
                                 quiet=TRUE)

    print("Getting candidates from FootprintFilter")
    total.time <- (system.time(tbl.out <- getCandidates(fp.filter)))[[3]]

    return(total.time)
}

timeRPsql <- function(host){

    print("Connecting through RPsql")
    conn.time <- system.time(db <- dbConnect(PostgreSQL(),
                                             user = "trena",
                                             password = "trena",
                                             dbname = "brain_hint_20",
                                             host = host))[[3]]
    
    # Make the 2 queries
    query.p0 <- "select loc, chrom, start, endpos from regions"
    query.p1 <- sprintf("where chrom='%s' and start >= %d and endpos <= %d", chrom, start, end)
    query.regions <- paste(query.p0, query.p1)
    
    reg.time <- system.time(tbl.reg <- DBI::dbGetQuery(db, query.regions))[[3]]
    
    loc.set <- sprintf("('%s')", paste(tbl.reg$loc, collapse="','"))
    query.hits <- sprintf("select * from hits where loc in %s", loc.set)
    
    hit.time <- system.time(tbl.hits <- DBI::dbGetQuery(db, query.hits))[[3]]
    
    dbDisconnect(db)

    return(list(conn = conn.time, reg = reg.time, hit = hit.time))
}


# Run both functions for each host 10 times and store in a data frame

# Create 10-rep sets of the 2 hosts
bdds <- rep("bddsrds.globusgenomics.org", 10)
khal <- rep("khaleesi", 10)

# Run the FpFilter method 10 times for bdds
fp.bdds <- sapply(bdds, timeFpFilter)

# Run the FpFilter method 10 times for khaleesi
fp.khal <- sapply(khal, timeFpFilter)

# Run the RPsql method 10 times for bdds and khaleesi
ps.bdds <- lapply(bdds, timeRPsql)
ps.khal <- lapply(khal, timeRPsql)

# More here to process the databases
