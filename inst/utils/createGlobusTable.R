library(htmlTable)
library(dplyr)

# Function for grabbing the list of databases currently on BDDS
getAvailableDBs <- function(host){
  
  db <- dbConnect(PostgreSQL(), 
                  user = "trena", 
                  password = "trena", 
                  host = host, 
                  dbname = "hg38")
  
  existing.databases <- dbGetQuery(db, "select datname from pg_database")[,1]
  
  # Pull out the databases I want using a grep
  important.dbs <- grep("(hint|wellington)_(20|16)",existing.databases, value = TRUE)
  dbDisconnect(db)
  return(important.dbs)
}

findApproxHits <- function(dbname, host){
  
  # Specify the connection
  db <- dbConnect(PostgreSQL(),
                  user = 'trena',
                  password = 'trena',
                  dbname = dbname,
                  host = host
  )
  # Get approximate hits
  numHits <- dbGetQuery(db, "select reltuples from pg_class where relname = 'hits'")[1,1]
  dbDisconnect(db)
  return(numHits)
}

# Run the 2 functions to get the table info
all.dbs <- getAvailableDBs(host = "bddsrds.globusgenomics.org")
all.hits <- sapply(all.dbs, findApproxHits, host = "bddsrds.globusgenomics.org")

# Assemble a data frame
df <- data_frame(Name = all.dbs,
                 Hits.in.millions = formatC(all.hits/1e6, digits = 1, format = "f"),
                 Tissue = gsub("(.*)_(hint|wellington).*","\\1",Name),
                 Method = gsub(".*_(hint|wellington).*","\\1",Name),
                 Seed = gsub(".*_(16|20)$", "\\1", Name))

df <- arrange(df, Tissue, Method, Hits.in.millions)

# https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html

html.txt <- htmlTable(df)

writeLines(html.txt, "./currentDBs.html")