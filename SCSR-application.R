library(glue)
library(devtools)

baseDir  <- '/data/villemin'
package <- '/data2/villemin/SingleCellSignalR/SingleCellSignalR'

devtools::install(glue("{baseDir}{package}"))
#library(SingleCellSignalR)
suppressPackageStartupMessages(library(SingleCellSignalR))

print("CreateDatabase on Request")
#createDatabase()

print("checkLastVersion FALSE")
#checkLastVersion(update=FALSE)

print("checkLastVersion TRUE")
#checkLastVersion(update=TRUE)

print("getInteractions")

interactions <- getInteractions()
head(interactions)
dim(interactions)

getInteractions(idRelease=1)
head(interactions)
dim(interactions)


getInteractions(idRelease=65)
