#' @title Interaction Dot Plot
#' @description Displays the level of expression of a LR pair in each cell.
#'
#' @details
#' `interaction` must be a number that corresponds to an interaction in the
#' SCSRInference object.
#' @details
#' if `pathway` is set then only receptors involved in those pathways will
#' be plotted.
#' @details
#' if `pathway` is set then the `dm` argument needs to be specified to give
#' information over the dataset.
#'
#' @param obj an SCSRInference object
#' @param dm an SCSRDataModel object
#' @param interaction the identifier of the interaction number in the 
#' interaction dataframe.
#' @param pathway the name of a pathway of interest.
#'
#' @return The function returns a R plot.
#' @import ggplot2
#' @export
#'
#' @examples
#' message('interactionDotPlot')
#' message('--dataPrepare')
#' data <- matrix(runif(1000,0,1),nrow=50,ncol=20)
#' rownames(data) <- paste('gene',seq_len(50))
#' obj <- dataPrepare(data)
#' message('--add Clustering')
#' obj <- addClustering(obj,cluster.id = sample(1:5,size = 20,replace = TRUE))
#' obj.int <- cellSignaling(obj)
#' #interactionDotPlot(obj.int, interaction = 1)

interactionDotPlot <- function(obj, dm = NULL, interaction, pathway = NULL) {


    if (!is(obj, "SCSRInference")) {
        stop("obj must be a SCSRInference object")
    }

    signal <- obj@LRinter

    if (interaction > length(signal)|length(interaction)>1){
        stop("interaction must refer to the indice of an interaction in the
            SCSRInference object.")
    }


    signal <- signal[[interaction]]
    name <- paste(names(signal[1])," --> ",names(signal)[2])
    names(signal) <- c("lig", "rec", "Type", "score")

    if (!is.null(pathway)){
        if (is.null(dm)) {
            stop("dm must be set when specifying a pathway")
        }else if (!is(dm, "SCSRDataModel")) {
            stop("obj must be a SCSRDataModel object")
        }

        data(PwC_ReactomeKEGG)
        if (!any(grepl(paste0(pathway, collapse = "|"), 
            PwC_ReactomeKEGG$pathway))) {
            stop(paste0(pathway, collapse = ","), "doesn't  
                correspond to any pathway in database. Please check your   
                spelling, use function findPathway().")
        } else {
            genes.to.plot <- unique(c(PwC_ReactomeKEGG[grepl(paste0(pathway,
                collapse = "|"), PwC_ReactomeKEGG$pathway), 1], 
                    PwC_ReactomeKEGG[grepl(paste0(pathway,
                collapse = "|"), PwC_ReactomeKEGG$pathway), 2]))
            if (dm@initial.organism != "hsapiens") {
                ortho <- data.frame(Hsapiens = rownames(dm@ncounts$matrix), 
                    species = dm@ncounts$initial.orthologs)
                genes.to.plot <- data.frame(Hsapiens = genes.to.plot)
                genes.to.plot <- merge(genes.to.plot, ortho, by.x = "Hsapiens",
                    order = FALSE)
                genes.to.plot <- genes.to.plot$species
            }
            signal <- signal[signal$rec%in%genes.to.plot,]
            if (nrow(signal)==0) message("No receptors involved in ",
                paste0(pathway, collapse = ","), " in the interaction ", name)

        }
    }

    p <- ""

    if (nrow(signal)>0){
        p <- ggplot(signal, aes(x=lig, y = rec, color = score, size = score)) +
                geom_point() + 
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                    hjust=1)) +
                ggtitle(name) +
                theme(axis.ticks = element_blank()) +
                scale_color_gradient(low="grey", high="darkblue")
    }

    p
}
