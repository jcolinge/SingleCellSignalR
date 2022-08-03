#' @title find pathways
#' @description Find pathway names within the database to use in functions.
#'
#' @details The `keyword` parameter is set to find pathways that contain it.
#'
#' @param keyword a character string
#'
#' @return The function returns a dataframe containing pathway names that 
#' contain the specified keyword.
#'
#' @export
#'
#' @examples
#' message('findPathways')
#' df <- findPathways('immune')

findPathways <- function(keyword) {

    data(PwC_ReactomeKEGG)

    pw.names <- unique(unlist(strsplit(PwC_ReactomeKEGG$pathway, ";")))

    pw.names <- pw.names[grepl(keyword, pw.names, ignore.case = TRUE)]

    message(length(pw.names), " pathways containing the keyword ",
        keyword, ".\n")

    pw.names <- data.frame(pathways = pw.names)

    message(paste0(capture.output(head(pw.names)$pathways), collapse = "\n"))

    return(pw.names)

}


