library(RMySQL)
library(RMariaDB)

packageVersion("RMySQL")
packageVersion("RMariaDB")

driver.1 <- RMySQL::MySQL()
driver.2 <- RMariaDB::MariaDB()

host <- "genome-mysql.cse.ucsc.edu"
user <- "genome"
dbname <- "hg38"

db.1 <- DBI::dbConnect(driver.1, user = user, host = host, dbname = dbname)
all.tableNames.1 <- DBI::dbListTables(db.1);
head(all.tableNames.1)
db.2 <- DBI::dbConnect(driver.2, user = user, host = host, dbname = dbname)
all.tableNames.2 <- DBI::dbListTables(db.2)
head(all.tableNames.2)

