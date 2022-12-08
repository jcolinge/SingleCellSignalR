library(glue)
library(devtools)

baseDir  <- '/data/villemin'
package <- '/data2/villemin/SingleCellSignalR/SingleCellSignalR'

devtools::install(glue("{baseDir}{package}"))

suppressPackageStartupMessages(library(SingleCellSignalR))

print("CreateDatabase on Request")
database::createDatabase()

print("checkLastVersion FALSE")
database::checkLastVersion(update=FALSE)

print("checkLastVersion TRUE")
database::checkLastVersion(update=TRUE)





