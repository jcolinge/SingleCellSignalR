#' Fetch the database from internet.
#'
#' Fetch LR database from remote location.
#'
#' @param onRequest logical True if you force
#' download again. This will overwrite 
#' pre-existing database. Default is True.
#'
#' @import httr
#' @importFrom cli col_cyan
#' @export
#' @examples
#' print("createDatabase")
#' createDatabase()
createDatabase <- function(onRequest=TRUE){

    #Default directory
    cacheDir <- Sys.getenv("SingleCellSignalR_CACHEDIR")
 
    if (!dir.exists(cacheDir))  
        dir.create(cacheDir)        
       
    url <-  Sys.getenv("SingleCellSignalR_DB_URL")
 
    databaseFilePath <- paste(cacheDir
        ,basename(url)
        ,sep = "/")

    if(!file.exists(databaseFilePath) | onRequest) {

        isDownloaded <- .downloadDatabase(url,databaseFilePath)
        if(!isDownloaded)
            stop("Ligand-Receptor database was not downloaded successfully.")
    else {
         cat(cli::cli_alert_info("Ligand-Receptor database downloaded with success.","\n"))
        }

    }

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

    databaseFilePathCopy <- paste(Sys.getenv("SingleCellSignalR_CACHEDIR")
        ,"SingleCellSignalR.copy.db"
        ,sep = "/")

    isDownloaded <- .downloadDatabase(Sys.getenv("SingleCellSignalR_DB_URL"),databaseFilePathCopy)
    if(!isDownloaded)
            stop("New Ligand-Receptor database was not downloaded successfully.")

    connexionObject <- DBI::dbCanConnect(RSQLite::SQLite(), databaseFilePathCopy)

    .checkDatabaseValidity(connexionObject)
    
    md5Local <- tools::md5sum(databaseFilePathCopy)

    SingleCellSignalRCon <- DBI::dbConnect(RSQLite::SQLite(), databaseFilePathCopy)
 
    md5DB <- DBI::dbGetQuery(SingleCellSignalRCon, "SELECT md5 FROM Release ORDER BY id DESC LIMIT 1")

    # If DB file is not the correct last database :
    if(md5Local!=md5DB){
         cat(cli::cli_alert_info("A new Ligand-Receptor database version is available !","\n") )
         if (update){ 
            #createDatabase(onRequest=update) 
            databaseFilePath <- paste(Sys.getenv("SingleCellSignalR_CACHEDIR")
            ,basename(Sys.getenv("SingleCellSignalR_DB_URL"))
            ,sep = "/")
            unlink(databaseFilePath)
            file.rename(databaseFilePathCopy,databaseFilePath)
            cat(cli::cli_alert("Ligand-Receptor database has been be updated.","\n"))

         }
         else { cat(cli::cli_alert("Install new version using `createDatabase(onRequest=TRUE)`.","\n"))} 
    }

    DBI::dbDisconnect(SingleCellSignalRCon)

}
     