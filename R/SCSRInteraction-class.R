library(methods)

#' SingleCellSignalR Interaction Object
#'
#' An S4 class to represent ligand-receptor interactions.
#'
#' @slot LRinter  A list of data frames describing the (ligand,receptor) doubles
#' with interaction type and LR score for each subpopulation combination.
#' @slot LRsubpop  A data frame describing the number of interactions
#' for each subpopulation combination.
#' @slot ligands   A vector of ligands, one entry per LR interaction.
#' @slot receptors   A vector of receptors, one entry per LR interaction.
#' @slot param   A list containing the calculation details of interaction.
#'
#' @details This class is a container for LR interactions along with
#' their statistical confidence. 
#' @export
#' @examples
#' new("SCSRInteraction")
#'
setClass("SCSRInteraction",
         slots=c(LRinter = "list",
                 LRsubpop = "data.frame",
                 ligands = "character",
                 receptors = "character",
                 param = "list"),
         prototype=list(
             LRinter=list(AB = data.frame(L = "x", R = "y", int.type = "paracrine",
                                LRscore = 0.5, stringsAsFactors = FALSE),
                          CD = data.frame(L = "v", R = "w", int.type = "paracrine",
                                LRscore = 0.5, stringsAsFactors = FALSE)),
             LRsubpop = data.frame(Subpop = c("AB","CD"), Number = c(1,1)),
             ligands = c("x","v"),
             receptors = c("y","w"),
             param = list(s.score = 0.5, logFC = 2, tol = 0.1)
         ))

setValidity("SCSRInteraction",
            function(object) {
                if(!is.list(object@LRinter))
                    return("LRinter is not a list")
                if(!is.data.frame(object@LRsubpop))
                    return("LRsubpop is not a dataframe")
                if(!is.character(object@ligands))
                    return("ligands is not a character")
                if(sum(unlist(lapply(object@LRinter,nrow)))!=length(object@ligands))
                    return("the length of ligands must be the sum of all interactions")
                if(!is.character(object@receptors))
                    return("receptors is not a character")
                if(length(object@receptors)!=length(object@ligands))
                    return("receptors must be the same length as ligands")
                if(!is.list(object@param))
                    return("param is not a list")

                TRUE
            }
)

setMethod("show", "SCSRInteraction",
          function(object) {
              cat("Interaction type: ", object@LRinter[[1]]$int.type[1], "\n", sep="")
              cat("Parameters: \n- S.score: ", param@s.score,"\n- logFC: ", 
                param@logFC, "\n- tolerance: ", param@tol, "\n", sep="")
              cat("LR interaction: \n")
              print(head(object@LRinter[order(object@LRinter$LRscore),
                        c("L", "R", "LRscore")]))
              cat("Subpopulation interaction counts:\n")
              print(object@LRsubpop)
          }
)


# Accessors & setters ========================================================

if (!isGeneric("LRinter")) {
    if (is.function("LRinter"))
        fun <- LRinter
    else
        fun <- function(x) standardGeneric("LRinter")
    setGeneric("LRinter", fun)
}
#' LRinter accessor
#' @export
setMethod("LRinter", "SCSRInteraction", function(x) x@LRinter)

if (!isGeneric("LRsubpop")) {
    if (is.function("LRsubpop"))
        fun <- LRsubpop
    else
        fun <- function(x) standardGeneric("LRsubpop")
    setGeneric("LRsubpop", fun)
}
#' LRsubpop accessor
#' @export
setMethod("LRsubpop", "SCSRInteraction", function(x) x@LRsubpop)

if (!isGeneric("ligands")) {
    if (is.function("ligands"))
        fun <- ligands
    else
        fun <- function(x) standardGeneric("ligands")
    setGeneric("ligands", fun)
}
#' ligands accessor
#' @export
setMethod("ligands", "SCSRInteraction", function(x) x@ligands)

if (!isGeneric("receptors")) {
    if (is.function("receptors"))
        fun <- receptors
    else
        fun <- function(x) standardGeneric("receptors")
    setGeneric("receptors", fun)
}
#' receptors accessor
#' @export
setMethod("receptors", "SCSRInteraction", function(x) x@receptors)


if (!isGeneric("param")) {
    if (is.function("param"))
        fun <- param
    else
        fun <- function(x) standardGeneric("param")
    setGeneric("param", fun)
}
#' Parameters accessor
#' @export
setMethod("param", "SCSRInteraction", function(x) x@param)
