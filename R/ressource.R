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
           

#' Import pathways from a file to dataframe 
#'
#' Pathways are defined in Reactome and
#' GoBP databases.
#' Those can be updated using
#' json files from
#' the Human Molecular Signatures Database (MSigDB)
#' at \url{https://www.gsea-msigdb.org/}
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
   
        if (! file.exists(file))
           stop("This file doesn't exist.")
        
        if(fileType=="json")
           db <- .formatPathwaysFromJson(file=file,
                pathwaySource=pathwaySource)
        else if(fileType=="gmt")
           db <- .formatPathwaysFromGmt(file=file,
                pathwaySource=pathwaySource)
        else {  stop("File format is not defined correctly.")}

        # Due to the fact React and Go are organized differently
        if(pathwaySource=="REACTOME"){
            names(db)<-c('Reactome name','Gene name','Reactome ID')
        }

        if(pathwaySource=="GO-BP")
            names(db)<-c('GO ID','Gene name','GO name')

return (db)

} #updatePathwaysFromFile

#' Format dataframe according to json input
#'
#' @param file    Path to file.
#' @param pathwaySource    Two options "GO-BP" or "REACTOME".
#'
#' @return Dataframe with pathwayID, geneName and pathwayName
#'
#' @importFrom foreach %do% %dopar%
#' @import doParallel
#' @import jsonlite
.formatPathwaysFromJson <- function(file,
    pathwaySource=NULL) {
    
    data <- jsonlite::read_json(file,simplifyVector = TRUE)

    db <-foreach::foreach(indexPathway=1:length(data), 
                    .combine = 'rbind') %dopar% { 
            data.frame(pathwayID=rep(data[[indexPathway]]$exactSource,length(data[[indexPathway]]$exactSource)),
                      geneName=unlist(data[[indexPathway]]$geneSymbols)[1:length(unlist(data[[indexPathway]]$geneSymbols))],
                      pathwayName= rep(names(data)[[indexPathway]],length(names(data)[[indexPathway]])))
    }

return (db)

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
#' \url{https://reactome.org/download/current/ReactomePathways.gmt.zip}
#' Note that you need to unzip the file to read the content.

#' The code is inspired from read.gmt function
#' from the gsa R package.
#'
#' @param file    Path to GMT file
#' @param pathwaySource   Two options "GO-BP" or "REACTOME"
#'
#' @return Dataframe with pathwayID, geneName and pathwayName
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


return(dataframeFromGmt)

} #.formatPathwaysFromGmt

 


#' Creache all ressources.
#'
#' Create cache for all ressources (pathways, or PWC network)
#' downloaded from the web when library is first loaded.
#' This part is handled with BiocFileCache.
#' Otherwise datatabase, is handled by another process
#' not relying on BiocFileCache instance.
#'
#' @param onRequest logical True if you force
#' download again. This will overwrite 
#' pre-existing database. Default is True.
#' @param verbose Default is FALSE
#'
#' @export 
#' @examples
#' if(FALSE)
#'  createRessources() 
createRessources <- function(onRequest=TRUE,verbose=FALSE) {

   cacheDir <-     Sys.getenv("SingleCellSignalR_CACHEDIR")
   ressourcesCacheDir <- paste(cacheDir,"ressources",sep="/")
     
   # Do it once, onLoad
   if(!dir.exists(ressourcesCacheDir) | onRequest) {
        .addCache(url=Sys.getenv("SingleCellSignalR_GO_URL"),cacheDir=ressourcesCacheDir,ressourceName="GO-BP",verbose=verbose)
        .addCache(url=Sys.getenv("SingleCellSignalR_Reactome_URL"),cacheDir=ressourcesCacheDir,ressourceName="Reactome",verbose=verbose)
        .addCache(url=Sys.getenv("SingleCellSignalR_PwC_URL"),cacheDir=ressourcesCacheDir,ressourceName="PwC",verbose=verbose)

    }
}


#' Add cache for ressources.
#'
#' Add cache for ressources (pathways, or PWC network)
#' downloaded from the web.
#' This part is handled with BiocFileCache.
#' Otherwise datatabase, is handled by another process
#' not relying on BiocFileCache instance.
#'
#' @param url    Path to file on the web.
#' @param ressourceName   Ressource name.
#' @param cacheDir   Absolute path to cache directory.
#' @param verbose   Default FALSE

#' @import BiocFileCache
#' @import httr
#' @keywords internal
.addCache <- function(url,cacheDir,ressourceName,verbose=FALSE) {

        bfc <- BiocFileCache::BiocFileCache(cacheDir,ask = FALSE)
   
        config <- httr::set_config(config(ssl_verifypeer = 0L,ssl_verifyhost = 0L))
            
        # if fname="exact" remove the unique identifier
        BiocFileCache::bfcadd(bfc,rname=ressourceName,config=config,fpath=url)

        if(verbose){
            print(BiocFileCache::bfccache(bfc))
            print(length(bfc)) 
            print(BiocFileCache::bfcinfo(bfc))
        }
}

#' Get ressource from the cache.
#'
#' Get  ressources (pathways, or PathwayCommons network 
#' from \url{https://www.pathwaycommons.org/})
#' stored in the cache.
#'
#' @param ressourceName   Ressource name.
#'
#' @importFrom cli cli_alert_danger
#' @export
#' @examples
#' reactome <-  getRessource(ressource="Reactome")
getRessource <- function(ressourceName=NULL) {

    if (!ressourceName %in% c("GO-BP","Reactome","PwC")){
        cat(cli::cli_alert_danger("GO-BP and Reactome are the only keywords alllowed.","\n"))
        stop() 
    }
    cacheDir <-     Sys.getenv("SingleCellSignalR_CACHEDIR")
    ressourcesCacheDir <- paste(cacheDir,"ressources",sep="/")

    bfc <- BiocFileCache::BiocFileCache(ressourcesCacheDir,ask = FALSE)

    dataframe <- .readFromCache(bfc=bfc,ressourceName=ressourceName)

    return(dataframe)
}

#' Read from the cache.
#'
#' Access  ressources (pathways, or PathwayCommons network 
#' from \url{https://www.pathwaycommons.org/})
#' stored in the cache.
#'
#' @param bfc Object of class BiocFileCache, created by a call to 
#' BiocFileCache::BiocFileCache()
#' @param ressourceName keyword associated to a specific ressourceName
#' @keywords internal
.readFromCache <- function(bfc,ressourceName) {

    cacheHits <- bfcquery(bfc,query=ressourceName,field="rname")
    if(nrow(cacheHits) == 0) {
        cat(cli::cli_alert_danger("No cache result found.","\n"))
        stop() 
    }
    else if(nrow(cacheHits) > 1) {
         cat(cli::cli_alert_danger("Multiple cache results found.","\n"))
         stop("Please clear your cache by running cacheClear()!")
    } else {
        rid <- cacheHits$rid
        result <- readRDS( bfc[[ rid ]] )
        return(result)
    }
}

#' Check existence of a record in the cache.
#'
#' Check if the cache record exists or not, by passing
#' to the function an associated keyword
#' associated to the ressource we are looking for.
#'
#' @param bfc Object of class BiocFileCache, created by a call to 
#' BiocFileCache::BiocFileCache()
#' @param ressourceName keyword associated to a specific ressource name
#' 
#' @keywords internal
#' @return logical This function returns TRUE if a record with 
#' the requested keyword already  exists in the file cache,
#'  otherwise returns FALSE.
.checkInCache <- function(bfc,ressourceName) {
    cacheHits <- bfcquery(bfc, query = ressourceName, field = "rname")
    as.logical(nrow(cacheHits))
}

#' Check valid cache.
#'
#' This function checks if a cache entry is a valid RDS file.
#' Returns TRUE if the cache entry is valid, FALSE otherwise.
#' In the case of an invalid file the cache entry and file are 
#' deleted.
#'
#' @param bfc Object of class BiocFileCache, created by a call to 
#' BiocFileCache::BiocFileCache()
#' @param ressourceName keyword associated to a specific ressource name.
#' @importFrom cli cli_alert_danger
#' @importFrom BiocFileCache bfcremove
#' @keywords internal
.checkValidCache <- function(bfc, ressourceName) {
    cacheHits <- bfcquery(bfc,query=ressourceName,field="rname")
    if(nrow(cacheHits) == 0) {
       cat(cli::cli_alert_danger("No cache result found.","\n"))
       stop() 
    }
    else if(nrow(cacheHits) > 1) {
         cat(cli::cli_alert_danger("Multiple cache results found.","\n"))
         stop("Please clear your cache by running cacheClear()!")
    } else {
        test <- tryCatch(is.list(infoRDS(cacheHits$rpath[1])), 
                         error = function(e) { return(FALSE) })
        if(!test) 
            BiocFileCache::bfcremove(bfc, cacheHits$rid[1])
        return(test)
    }
}


#' Delete cache content.
#'
#' Delete the content of cache directory.
#'
#' @importFrom BiocFileCache removebfc
#' @export
#' @examples
#  if(FALSE)
#   cacheClear()
cacheClear <- function() {

    cacheDir <-     Sys.getenv("SingleCellSignalR_CACHEDIR")
    ressourcesCacheDir <- paste(cacheDir,"ressources",sep="/")

    bfc <- BiocFileCache::BiocFileCache(ressourcesCacheDir, ask = FALSE)
    BiocFileCache::removebfc(bfc, ask = FALSE)

    #dir.create(ressourcesCacheDir)
    message("SingleCellSignalR cache has been deleted.\n", 
                "- Location: ", ressourcesCacheDir, "\n",
                "- No. of files: 0", "\n")

}

#' Get cache content informations..
#'
#' Get cache content informations.
#'
#' @importFrom BiocFileCache
#' @importFrom cli cli_alert_danger
#' @export
#' @examples
#  if(FALSE)
#   cacheInfo()
cacheInfo <- function() {

    cacheDir <-  Sys.getenv(x = "SingleCellSignalR_CACHEDIR")
    ressourcesCacheDir <- paste(cacheDir,"ressources",sep="/")
    
    # safeguard
    if(!dir.exists(ressourcesCacheDir)) {
        dir.create(ressourcesCacheDir)
    } 

    files <-  list.files(ressourcesCacheDir)

    if(length(files)==0) {
        message("SingleCellSignalR cache uninitialized.\n", 
                "- Location: ", ressourcesCacheDir, "\n",
                "- No. of files: ", length(files), "\n")

    } else {
        
        bfc <- BiocFileCache::BiocFileCache(ressourcesCacheDir, ask = FALSE)
        files <- bfcinfo(bfc)$rpath
        total_size <- sum(file.size(files))
        size_obj <- structure(total_size, class = "object_size")
    
        message("SingleCellSignalR cache: \n", 
                "- Location: ", ressourcesCacheDir, "\n",
                "- No. of files: ", length(files), "\n",
                "- Total size: ", format(size_obj, units = "auto"), "\n")
    }

    return(invisible(ressourcesCacheDir))
}

