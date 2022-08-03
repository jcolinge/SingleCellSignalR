#' Ligand/Receptor interactions data table
#'
#' @format A data frame with 3251 rows of 13 variables:
#' \describe{
#'     \item{ligand}{ligand gene symbol}
#'     \item{receptor}{receptor gene symbol}
#'     \item{source}{provenance of the interaction}
#'     \item{PMIDs}{PubmedID of the publication reporting the interaction}
#'     ...
#' }
#' @source \url{EDF R\&D}
#' @usage data(LRdb)
#' @docType data
#' @keywords internal, dataset
"LRdb"

#' Pathway Commons Reactome KEGG 2019-05-08
#'
#' @format A data frame with 26067 rows of 6 variables:
#' \describe{
#'     \item{a.gn}{interactant 1 gene symbol}
#'     \item{b.gn}{interactant 2 gene symbol}
#'     \item{type}{interaction type}
#'     \item{pathway}{Associated pathway}
#'     ...
#' }
#' @source \url{EDF R\&D}
#' @usage data(PwC_ReactomeKEGG)
#' @docType data
#' @keywords internal, dataset
"PwC_ReactomeKEGG"

#' A list of cell types markers
#'
#' @format A data frame with 95 rows of 15 variables:
#' \describe{
#' \item{TNBC}{Triple Negative Breast Cancer markers}
#'     \item{HER+}{HER+ Breast Cancer markers}
#'     \item{ER+}{ER+ Breast Cancer markers}
#'     \item{T-cells}{T cells markers}
#'     \item{B-cells}{B cells markers}
#'     \item{Macrophages}{Macrophages markers}
#'     \item{Endothelial cells}{Endothelial cells markers}
#'     \item{CAFs}{Cancer Associated Fibroblasts markers}
#'     \item{melanoma}{Melanoma markers}
#'     \item{Cytotoxic cells}{Cytotoxic cells markers}
#'     \item{DC}{Dendritic cells markers}
#'     \item{Mast cells}{Mastocytes cells markers}
#'     \item{Neutrophils}{Neutrophils cells markers}
#'     \item{NK cells}{Natural Killer cells markers}
#'     \item{Treg}{Regulatory T-cells markers}
#'     ...
#' }
#' @source \url{EDF R\&D}
#' @usage data(markers_default)
#' @docType data
#' @keywords internal, dataset
"markers_default"

#' Example dataset
#'
#' @format A data frame with 1520 rows of 401 variables:
#'
#' @source \url{EDF R\&D}
#' @usage data(example_dataset)
#' @docType data
#' @keywords dataset
"example_dataset"

#' Example dataset - Patient80
#'
#' @format A data frame with 500 rows of 480 variables:
#'
#' @source \url{EDF R\&D}
#' @usage data(Patient80)
#' @docType data
#' @keywords dataset
"Patient80"

#' Example dataset - Mouse
#'
#' @format A data frame with 500 rows of 425 variables:
#'
#' @source \url{EDF R\&D}
#' @usage data(Mouse)
#' @docType data
#' @keywords dataset
"Mouse"