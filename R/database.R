#' Fetch the database from internet.
#'
#' Fetch LR database from remote location.
#'
#' @param onRequest logical True if you force
#' download again. This will overwrite 
#' pre-existing database. Default is True.
#' @param verbose Logical TRUE/FALSE
#'
#' @import httr
#' @importFrom cli col_cyan
#' @export
#' @examples
#' print("createDatabase")
#' createDatabase()
createDatabase <- function(onRequest=TRUE,verbose=FALSE){

    #Default directory 
    cacheDir <- Sys.getenv("SingleCellSignalR_CACHEDIR")
    databaseCacheDir <- paste(cacheDir,"database",sep="/")
    
    if(!dir.exists(databaseCacheDir))
        dir.create(databaseCacheDir,recursive = TRUE)

    url <-  Sys.getenv("SingleCellSignalR_DB_URL")
 
    databaseFilePath <- paste(databaseCacheDir
        ,basename(url)
        ,sep = "/")

    if(!file.exists(databaseFilePath) | onRequest) {

        isDownloaded <- .downloadDatabase(url,databaseFilePath)
        if(!isDownloaded)
            stop("Ligand-Receptor database was not downloaded successfully.")
        .addCache(fpath=databaseFilePath,cacheDir=databaseCacheDir,resourceName="LRdb",verbose=verbose)
    }
    
    else {cli::cli_alert_info("{.val LRdb} database downloaded with success.","\n")}

    if(file.exists(databaseFilePath)) {
         
        connexionObject <- DBI::dbCanConnect(RSQLite::SQLite(), databaseFilePath)

        .checkDatabaseValidity(connexionObject=connexionObject)

    }

    return(invisible())
}


#' Fetch the database from internet.
#'
#' Fetch LR database from remote location.
#'
#' @param url File URL. 
#' @param databaseFilePath Path to database file. 
#'
#' @import httr rappdirs
#' @importFrom cli col_cyan
#' @NoRd
.downloadDatabase <- function(url,databaseFilePath){

    cat(cli::col_cyan("Download Ligand-Receptor database...","\n"))

    isValid <- TRUE
        httr::set_config(config(ssl_verifypeer = 0L,ssl_verifyhost = 0L))
        response <- httr::GET(url, 
            httr::write_disk(databaseFilePath, overwrite=TRUE),
            httr::progress())

        if (httr::http_error(response)) {
            unlink(databaseFilePath)
            httr::warn_for_status(response, paste("find data at", url))
            isValid <- FALSE
        }
   return(isValid)
}

#' Check validity of database
#'
#' Control connexion is ok.
#'
#' @param connexionObject DBI::dbCanConnect object
#' to test if it's a valid connexion is possible.
#'
#' @import DBI RSQLite 
#' @importFrom cli cli_alert_danger
#' @NoRd
.checkDatabaseValidity <- function(connexionObject) {

    # check file is a database
    tryCatch(connexionObject,
            warning = function(e) {
              cat(cli::cli_alert_danger("Ligand-Receptor database is corrupted.","\n"))
              stop("The file provided is not a valid database.")
            }
          )
}

#' Check database version.
#'
#' Control if installed database is the last version.
#'
#' @param update logical Decide whether or not
#' to update the database. Default False.
#'
#' md5sum of database file is stored inside
#' the database.
#  If the localy recomputed md5sum string
#' is different, that means a new version
#' is available.
#' If update set to True, database is directly updated.
#' by calling \code{\link{createDatabase(onRequest=TRUE)}}
#'
#' @import httr DBI RSQLite tools
#' @importFrom cli cli_alert_danger cli_alert_info cli_alert
#' @export
#' @examples
#' print("checkLastVersion")
#' checkLastVersion(update=FALSE)
checkDatabaseLastVersion <- function(update=FALSE) {

 

}
     