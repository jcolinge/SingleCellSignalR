#' Retrieve LR complexes
#'
#' Fetch LR complexes from database and
#' and return a dataframe
#'
#' @param idRelease integer id version Release
#' Default is NULL so last version is selected.
#'
#' @import DBI RSQLite
#' @export
#' @examples
#' print("getComplexes")
#' getComplexes()
getComplexes <- function(idRelease=NULL) { 

    print("getComplexes")

    cacheDir <- Sys.getenv("SingleCellSignalR_CACHEDIR")
    url      <-  Sys.getenv("SingleCellSignalR_URL")
 
    databaseFilePath <- paste(cacheDir,basename(url),sep = "/")

    SingleCellSignalRCon <- dbConnect(RSQLite::SQLite(), databaseFilePath)
    if(is.null(idRelease))
        release <- DBI::dbGetQuery(SingleCellSignalRCon, 'SELECT id FROM Release ORDER BY id DESC LIMIT 1')
 
    else {
        release <- DBI::dbGetQuery(SingleCellSignalRCon, 'SELECT id FROM Release WHERE id = ?',
        , params = list(idRelease))
    }
 
    if(nrow(release)==0)
        cat(cli::cli_abort("ID Release {idRelease} doesn't exist.","\n"))

    complexesL <- DBI::dbGetQuery(SingleCellSignalRCon, 
        'SELECT Clex.name, Clex.size,Clex.source,Comp.name,Comp.type
        FROM Complex as Clex  
        inner join Component as Comp
        inner join Complex_Component as CC
        on Comp.id = CC."id.component_fk" AND  Clex.id = CC."id.complex_fk" 
        where Comp."id.release_fk" = ? AND Comp.type = ?'
        , params = list(release$id,"L"))

    
    complexesR <- DBI::dbGetQuery(SingleCellSignalRCon, 
   'SELECT Clex.name, Clex.size,Clex.source,Comp.name,Comp.type
        FROM Complex as Clex  
        inner join Component as Comp
        inner join Complex_Component as CC
        on Comp.id = CC."id.component_fk" AND  Clex.id = CC."id.complex_fk" 
        where Comp."id.release_fk" = ? AND Comp.type = ?'
      , params = list(release$id,"R"))

    DBI::dbDisconnect(SingleCellSignalRCon)

    complexes <- rbind(complexesL, complexesR)

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
#' @import DBI RSQLite
#' @importFrom cli cli_abort
#' @export
#' @examples
#' print("getInteractions")
#' getInteractions()
getInteractions <- function(idRelease=NULL) { 

    cacheDir <- Sys.getenv("SingleCellSignalR_CACHEDIR")
    url      <-  Sys.getenv("SingleCellSignalR_URL")
 
    databaseFilePath <- paste(cacheDir,basename(url),sep = "/")

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
        'SELECT Comp.name, Comp.description,Inter.source 
        FROM Component as Comp inner join Interaction as Inter
         on Comp.id = Inter."id.ligand_fk" where Comp."id.release_fk" = ?'
        , params = list(release$id))

    colnames(ligands)[which(names(ligands)=="name")] <- "Ligand"               
    colnames(ligands)[which(names(ligands)=="description")] <- "Ligand.name"               
    colnames(ligands)[which(names(ligands)=="source")] <- "Ligand.source"               

    receptors <- DBI::dbGetQuery(SingleCellSignalRCon, 
    'SELECT Comp.name, Comp.description,Inter.source 
    FROM Component as Comp inner join Interaction as Inter
     on Comp.id = Inter."id.receptor_fk"  where Comp."id.release_fk" = ?'
        , params = list(release$id))

    colnames(receptors)[which(names(receptors)=="name")] <- "Receptor"               
    colnames(receptors)[which(names(receptors)=="description")] <- "Receptor.name"               
    colnames(receptors)[which(names(receptors)=="source")] <- "Receptor.source"     

    DBI::dbDisconnect(SingleCellSignalRCon)

    interactions <- cbind(ligands, receptors)

    return(invisible(interactions))
}
