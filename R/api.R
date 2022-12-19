#' Retrieve LR complexes
#'
#' Fetch LR complexes from database and
#' and return a dataframe
#'
#' @param idRelease integer id version Release
#' Default is NULL so last version is selected.
#'
#' @import DBI RSQLite BiocFileCache
#' @export
#' @examples
#' print("getComplexes")
#' getComplexes()
getComplexes <- function(idRelease=NULL) { 

    cacheDir <- Sys.getenv("SingleCellSignalR_CACHEDIR")
    databaseCacheDir <- paste(cacheDir,"database",sep="/")

    url      <-  Sys.getenv("SingleCellSignalR_DB_URL")
 
    bfc <- BiocFileCache::BiocFileCache(databaseCacheDir,ask = FALSE)
   
    databaseFilePath <- BiocFileCache::bfcrpath(bfc, rids = "BFC1")

    print(databaseFilePath)

    SingleCellSignalRCon <- dbConnect(RSQLite::SQLite(), databaseFilePath)

    if(is.null(idRelease))
        release <- DBI::dbGetQuery(SingleCellSignalRCon, 'SELECT id FROM Release ORDER BY id DESC LIMIT 1')
 
    else {
        release <- DBI::dbGetQuery(SingleCellSignalRCon, 'SELECT id FROM Release WHERE id = ?',
        , params = list(idRelease))
    }
 
    if(nrow(release)==0)
        cat(cli::cli_abort("ID Release {idRelease} doesn't exist.","\n"))

    complexes <- DBI::dbGetQuery(SingleCellSignalRCon, 
        'SELECT Clex.name,Clex.description,Clex.size,Comp.name,Comp.type,CC.stoichiometry,Clex.sources,Clex.pmids 
        FROM Complex as Clex  
        inner join Component as Comp
        inner join Complex_Component as CC
        on Comp.id = CC."id.component_fk" AND  Clex.id = CC."id.complex_fk" 
        where Comp."id.release_fk" = ?'
        , params = list(release$id))


    DBI::dbDisconnect(SingleCellSignalRCon)

    return(invisible(complexes))

}

#' Retrieve LR interactions.
#'
#' Fetch LR interactions from database and
#' and return a dataframe
#'
#' @param idRelease integer id version Release
#' Default is NULL so last version is selected.
#'
#' @import DBI RSQLite BiocFileCache
#' @importFrom cli cli_abort
#' @export
#' @examples
#' print("getInteractions")
#' getInteractions()
getInteractions <- function(idRelease=NULL) { 

    cacheDir <- Sys.getenv("SingleCellSignalR_CACHEDIR")
    url      <-  Sys.getenv("SingleCellSignalR_DB_URL")
 
    databaseCacheDir <- paste(cacheDir,"database",sep="/")
    
    bfc <- BiocFileCache::BiocFileCache(databaseCacheDir,ask = FALSE)
    databaseFilePath <- BiocFileCache::bfcrpath(bfc, rids = "BFC1")

    print(databaseFilePath)

    SingleCellSignalRCon <- dbConnect(RSQLite::SQLite(), databaseFilePath)
    if(is.null(idRelease))
        release <- DBI::dbGetQuery(SingleCellSignalRCon, 'SELECT id FROM Release ORDER BY id DESC LIMIT 1')
 
    else {
        release <- DBI::dbGetQuery(SingleCellSignalRCon, 'SELECT id FROM Release WHERE id = ?',
        , params = list(idRelease))
    }
 
    if(nrow(release)==0)
        cat(cli::cli_abort("ID Release {idRelease} doesn't exist.","\n"))
  
  ligands <- DBI::dbGetQuery(SingleCellSignalRCon, 
      'SELECT DISTINCT(Comp.id),Comp.name, Comp.description
      FROM Component as Comp inner join Interaction as Inter
       on Comp."id" = Inter."id.ligand_fk" where Comp."id.release_fk" = ?'
      , params = list(release$id))

  colnames(ligands)[which(names(ligands)=="id")] <- "id.ligand_fk"               
  colnames(ligands)[which(names(ligands)=="name")] <- "ligand"               
  colnames(ligands)[which(names(ligands)=="description")] <- "ligand.name"               

  receptors <- DBI::dbGetQuery(SingleCellSignalRCon, 
  'SELECT DISTINCT(Comp.id),Comp.name, Comp.description
  FROM Component as Comp inner join Interaction as Inter
   on Comp."id" = Inter."id.receptor_fk"  where Comp."id.release_fk" = ?'
      , params = list(release$id))

  colnames(receptors)[which(names(receptors)=="id")] <- "id.receptor_fk"               
  colnames(receptors)[which(names(receptors)=="name")] <- "receptor"               
  colnames(receptors)[which(names(receptors)=="description")] <- "receptor.name"    

  interactions <- DBI::dbGetQuery(SingleCellSignalRCon, 
  'SELECT Inter."id.ligand_fk", Inter."id.receptor_fk" ,Inter."sources",Inter."pmids"
  FROM Interaction as Inter inner join Component as Comp on Comp."id" = Inter."id.receptor_fk" where Comp."id.release_fk" = ?'
      , params = list(release$id))

  LRdb <- inner_join(receptors,interactions,  by='id.receptor_fk')
  LRdb <- inner_join(ligands , LRdb,  by='id.ligand_fk')
  LRdb$LR <- paste(LRdb$ligand,LRdb$receptor,sep="/")

  LRdb <- LRdb[,c("LR","ligand","ligand.name","receptor","receptor.name","sources","pmids")]
  
  DBI::dbDisconnect(SingleCellSignalRCon)

  return(invisible(LRdb))
}
