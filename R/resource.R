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
#' @param resourceName    Two options "GO-BP" or "REACTOME".
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
                        resourceName=NULL){

        fileType <- match.arg(fileType)

        if (! resourceName %in% c("GO-BP","REACTOME"))
           stop("GO-BP and REACTOME are the only keywords alllowed.")
   
        if (! file.exists(file))
           stop("This file doesn't exist.")
        
        if(fileType=="json")
           db <- .formatPathwaysFromJson(file=file,
                resourceName=resourceName)
        else if(fileType=="gmt")
           db <- .formatPathwaysFromGmt(file=file,
                resourceName=resourceName)
        else {  stop("File format is not defined correctly.")}

        # Due to the fact React and Go are organized differently
        if(resourceName=="REACTOME"){
            names(db)<-c('Reactome name','Gene name','Reactome ID')
        }

        if(resourceName=="GO-BP")
            names(db)<-c('GO ID','Gene name','GO name')

return (db)

} #updatePathwaysFromFile

#' Format dataframe according to json input
#'
#' @param file    Path to file.
#' @param resourceName    Two options "GO-BP" or "REACTOME".
#'
#' @return Dataframe with pathwayID, geneName and pathwayName
#'
#' @importFrom foreach %do% %dopar%
#' @import doParallel
#' @import jsonlite
.formatPathwaysFromJson <- function(file,
    resourceName=NULL) {
    
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
#' @param resourceName   Two options "GO-BP" or "REACTOME"
#'
#' @return Dataframe with pathwayID, geneName and pathwayName
#'
#' @importFrom foreach %do% %dopar%
#' @import doParallel
.formatPathwaysFromGmt<- function(file,
    resourceName=NULL){

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

 


#' Creache all resources.
#'
#' Create cache for all resources (pathways, or PWC network)
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
#'  createResources() 
createResources <- function(onRequest=TRUE,verbose=FALSE) {

   cacheDir <-     Sys.getenv("SingleCellSignalR_CACHEDIR")
   resourceName <- paste(cacheDir,"resources",sep="/")
     
   # Do it once, onLoad
   if(!dir.exists(resourcesCacheDir) | onRequest) {
        .addCache(url=Sys.getenv("SingleCellSignalR_GO_URL"),cacheDir=resourcesCacheDir,resourceName="GO-BP",verbose=verbose)
        .addCache(url=Sys.getenv("SingleCellSignalR_Reactome_URL"),cacheDir=resourcesCacheDir,resourceName="Reactome",verbose=verbose)
        .addCache(url=Sys.getenv("SingleCellSignalR_PwC_URL"),cacheDir=resourcesCacheDir,resourceName="PwC",verbose=verbose)

    }

return(invisible(NULL))

}


#' Add cache for resources.
#'
#' Add cache for resources (pathways, or PWC network)
#' downloaded from the web.
#' This part is handled with BiocFileCache.
#' Otherwise datatabase, is handled by another process
#' not relying on BiocFileCache instance.
#'
#' @param url    Path to file on the web.
#' @param resourceName   Ressource name.
#' @param cacheDir   Absolute path to cache directory.
#' @param verbose   Default FALSE

#' @import BiocFileCache
#' @import httr
#' @keywords internal
.addCache <- function(url,cacheDir,resourceName,verbose=FALSE) {

        bfc <- BiocFileCache::BiocFileCache(cacheDir,ask = FALSE)
   
        config <- httr::set_config(config(ssl_verifypeer = 0L,ssl_verifyhost = 0L))
            
        # if fname="exact" remove the unique identifier
        BiocFileCache::bfcadd(bfc,rname=resourceName,config=config,fpath=url)

        if(verbose){
            print(BiocFileCache::bfccache(bfc))
            print(length(bfc)) 
            print(BiocFileCache::bfcinfo(bfc))
        }

return(invisible(NULL))

}

#' Get ressource from the cache.
#'
#' Get  resources (pathways, or PathwayCommons network 
#' from \url{https://www.pathwaycommons.org/})
#' stored in the cache.
#'
#' @param resourceName   Ressource name.
#'
#' @importFrom cli cli_alert_danger
#' @export
#' @examples
#' reactome <-  getRessource(ressource="Reactome")
getRessource <- function(resourceName=NULL) {

    if (!resourceName %in% c("GO-BP","Reactome","PwC")){
        cat(cli::cli_alert_danger("GO-BP and Reactome are the only keywords alllowed.","\n"))
        stop() 
    }
    cacheDir <-     Sys.getenv("SingleCellSignalR_CACHEDIR")
    resourcesCacheDir <- paste(cacheDir,"resources",sep="/")

    bfc <- BiocFileCache::BiocFileCache(resourcesCacheDir,ask = FALSE)

    dataframe <- .readFromCache(bfc=bfc,resourceName=resourceName)

    return(dataframe)
}

#' Read from the cache.
#'
#' Access  resources (pathways, or PathwayCommons network 
#' from \url{https://www.pathwaycommons.org/})
#' stored in the cache.
#'
#' @param bfc Object of class BiocFileCache, created by a call to 
#' BiocFileCache::BiocFileCache()
#' @param resourceName keyword associated to a specific resourceName
#' @keywords internal
.readFromCache <- function(bfc,resourceName) {

    cacheHits <- bfcquery(bfc,query=resourceName,field="rname")
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
#' @param resourceName keyword associated to a specific ressource name
#' 
#' @keywords internal
#' @return logical This function returns TRUE if a record with 
#' the requested keyword already  exists in the file cache,
#'  otherwise returns FALSE.
.checkInCache <- function(bfc,resourceName) {
    cacheHits <- bfcquery(bfc, query = resourceName, field = "rname")
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
#' @param resourceName keyword associated to a specific ressource name.
#' @importFrom cli cli_alert_danger
#' @importFrom BiocFileCache bfcremove
#' @keywords internal
.checkValidCache <- function(bfc, resourceName) {
    cacheHits <- bfcquery(bfc,query=resourceName,field="rname")
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
    resourcesCacheDir <- paste(cacheDir,"resources",sep="/")

    bfc <- BiocFileCache::BiocFileCache(resourcesCacheDir, ask = FALSE)
    BiocFileCache::removebfc(bfc, ask = FALSE)

    #dir.create(resourcesCacheDir)
    message("SingleCellSignalR cache has been deleted.\n", 
                "- Location: ", resourcesCacheDir, "\n",
                "- No. of files: 0", "\n")

    return(invisible(NULL))

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
    resourcesCacheDir <- paste(cacheDir,"resources",sep="/")
    
    # safeguard
    if(!dir.exists(resourcesCacheDir)) {
        dir.create(resourcesCacheDir)
    } 

    files <-  list.files(resourcesCacheDir)

    if(length(files)==0) {
        message("SingleCellSignalR cache uninitialized.\n", 
                "- Location: ", resourcesCacheDir, "\n",
                "- No. of files: ", length(files), "\n")

    } else {
        
        bfc <- BiocFileCache::BiocFileCache(resourcesCacheDir, ask = FALSE)
        files <- bfcinfo(bfc)$rpath
        total_size <- sum(file.size(files))
        size_obj <- structure(total_size, class = "object_size")
    
        message("SingleCellSignalR cache: \n", 
                "- Location: ", resourcesCacheDir, "\n",
                "- No. of files: ", length(files), "\n",
                "- Total size: ", format(size_obj, units = "auto"), "\n")
    }

    return(invisible(resourcesCacheDir))
}

