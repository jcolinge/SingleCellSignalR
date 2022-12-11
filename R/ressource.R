#' Fetch the database from internet.
#'
#' Fetch LR database from remote location.
#'
#' @param onRequest logical True if you force
#' download again. This will overwrite 
#' pre-existing database. Default is True.
#'
#' @import httr rappdirs
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
       
    url <-  Sys.getenv("SingleCellSignalR_URL")
 
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

        .checkValidity(connexionObject=connexionObject)

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
.checkValidity <- function(connexionObject) {

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
checkLastVersion <- function(update=FALSE) {

    databaseFilePathCopy <- paste(Sys.getenv("SingleCellSignalR_CACHEDIR")
        ,"SingleCellSignalR.copy.db"
        ,sep = "/")

    isDownloaded <- .downloadDatabase(Sys.getenv("SingleCellSignalR_URL"),databaseFilePathCopy)
    if(!isDownloaded)
            stop("New Ligand-Receptor database was not downloaded successfully.")

    connexionObject <- DBI::dbCanConnect(RSQLite::SQLite(), databaseFilePathCopy)

    .checkValidity(connexionObject)
    
    md5Local <- tools::md5sum(databaseFilePathCopy)

    SingleCellSignalRCon <- DBI::dbConnect(RSQLite::SQLite(), databaseFilePathCopy)
 
    md5DB <- DBI::dbGetQuery(SingleCellSignalRCon, "SELECT md5 FROM Release ORDER BY id DESC LIMIT 1")

    # If DB file is not the correct last database :
    if(md5Local!=md5DB){
         cat(cli::cli_alert_info("A new Ligand-Receptor database version is available !","\n") )
         if (update){ 
            #createDatabase(onRequest=update) 
            databaseFilePath <- paste(Sys.getenv("SingleCellSignalR_CACHEDIR")
            ,basename(Sys.getenv("SingleCellSignalR_URL"))
            ,sep = "/")
            unlink(databaseFilePath)
            file.rename(databaseFilePathCopy,databaseFilePath)
            cat(cli::cli_alert("Ligand-Receptor database has been be updated.","\n"))

         }
         else { cat(cli::cli_alert("Install new version using `createDatabase(onRequest=TRUE)`.","\n"))} 
    }

    DBI::dbDisconnect(SingleCellSignalRCon)

}
           

#' Import pathways from a file to dataframe 
#'
#' Pathways are defined in Reactome and
#' GoBP databases.
#' Those can be updated using
#' json files from
#' the Human Molecular Signatures Database (MSigDB)
#' at \href{URL}{https://www.gsea-msigdb.org/}
#'
#' \code{resetDownstreamPathways} is a function
#' we provide to user to refresh REACTOME 
#' and GO-BP content included in BulkSignalR.
#'
#' Gmt file format also can be imported.
#'
#' @param file    Path to file.
#' @param fileType    Default is Json.
#' @param pathwaySource    Two options "GO-BP" or "REACTOME".
#'
#' @return dataframe with ID, Name and Gene
#'
#' @import jsonlite
#' @export
#' @examples
#' print('updatePathwaysFromFile')
#' if(FALSE)
#'    updatePathwaysFromFile(file,"GO-BP")
#'
updatePathwaysFromFile <- function(file,
                        fileType=c("json","gmt"),
                        pathwaySource=NULL){

        fileType <- match.arg(fileType)

        if (! pathwaySource %in% c("GO-BP","REACTOME"))
           stop("GO-BP and REACTOME are the only keywords alllowed.")
   
        if (! file.exists(jsonFile))
           stop("This file doesn't exist.")
        
        if(fileType=="json")
           db <- .formatPathwaysFromJson(file=jsonFile,
                pathwaySource=pathwaySource)
        else if(fileType=="gmt")
           db <- .formatPathwaysFromGmt(file=jsonFile,
                pathwaySource=pathwaySource)
        else {  stop("File format is not defined correctly.")}

return (db)

} #updatePathwaysFromFile

#' Format dataframe according to json input
#'
#' @param file    Path to file.
#' @param pathwaySource    Two options "GO-BP" or "REACTOME".
#'
.formatPathwaysFromJson <- function(file,
    pathwaySource=NULL) {
    
    d <- jsonlite::read_json(jsonFile,simplifyVector = TRUE)

    if(pathwaySource=="REACTOME")
       db <-  data.frame('Reactome ID'=character(),'Gene name'=character(),'Reactome name'=character()) 

    if(pathwaySource=="GO-BP")
       db <-  data.frame('GO ID'=character(),'Gene name'=character(),'GO name'=character()) 

} # .formatPathwaysFromJson

#' Transform gmt file to dataframe 
#'
#' We note discrepancy between format available 
#' over internet. 
#'
#' Here we consider a valid gmt file format defined
#' on each lines as follows :  
#' First is Pathway name, 
#' Then comes the ID,
#' Finally you will find genes symbols
#' according to the pathway defined on the line.
#'
#' You can find an example here.
#' - For Reactome. (Directly from their website)
#' \href{URL}{https://reactome.org/download/current/ReactomePathways.gmt.zip}
#' Note that you need to unzip the file to read the content.

#' The code is inspired from read.gmt function
#' from gsa R package.
#'
#' @param file    Path to GMT file
#' @param pathwaySource   Two options "GO-BP" or "REACTOME"
#'
#' @return gmt file as dataframe
#'
#' @importFrom foreach %do% %dopar%
#' @import doParallel
.formatPathwaysFromGmt<- function(file,
    pathwaySource=NULL){

        read1 <-scan(file,what=list("",""),sep="\t", 
            quote=NULL, fill=TRUE, flush=TRUE,multi.line=FALSE)

        read1 <- setNames(read1, c("ID","Description"))

        geneset.ids<-read1[1][[1]]
        geneset.descriptions<-read1[2][[1]]

        read2<-scan(file,what="",sep="\t", quote=NULL)
        read2<- read2[read2!=""]

        nn<-length(geneset.ids)
        n<-length(read2)
        ox<-rep(NA,nn)

        # Compute indice of pathway
        ii<-1
        for(i in 1:nn){
             while((read2[ii]!=geneset.descriptions[i]) | (read2[ii+1]!=geneset.ids[i]) ){
               ii=ii+1
            }
         ox[i]=ii   
         ii=ii+1
        }

        genesets=vector("list",nn)

        for(i in 1:(nn-1)){
          i1<-ox[i]+2
          i2<-ox[i+1]-1
          geneset.descriptions[i]<-read2[ox[i]+1]
          geneset.ids[i]<-read2[ox[i]] 
          genesets[[i]]<-read2[i1:i2]
        }

        geneset.ids[nn]<-read2[ox[nn]] 
        geneset.descriptions[nn]=read2[ox[nn]+1]
        genesets[[nn]]=read2[(ox[nn]+2):n]

        data=list(geneset.ids=geneset.ids,
                  geneset.descriptions=geneset.descriptions,
                  genesets=genesets
                 )

        dataframeFromGmt <-foreach::foreach(i= 1:length(data$geneset.ids), 
                    .combine = 'rbind') %dopar% { 

            data.frame(a=rep(data$geneset.ids[[i]],length(data$genesets[[i]])),
                      b=data$genesets[[i]][1:length(data$genesets[[i]])],
                      c=rep(data$geneset.descriptions[[i]],length(data$genesets[[i]])))
        }

        # Due to the fact React and Go are organized differently
        if(pathwaySource=="REACTOME"){
            names(dataframeFromGmt)<-c('Reactome name','Gene name','Reactome ID')
        }

        if(pathwaySource=="GO-BP")
            names(dataframeFromGmt)<-c('GO ID','Gene name','GO name')

return(dataframeFromGmt)

} #.formatPathwaysFromGmt
