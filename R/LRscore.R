#' Calculation of the LRscore
#'
#' @param l a value (or a vector) of the mean ligand expression in the
#' secreting cluster
#' @param r a value (or a vector) of the mean receptor expression in the
#' receiving cluster
#' @param s a value for scaling the score (usually the mean of the whole
#' read count table, the median or another similar value is possible),
#' must be over 0
#'
#' @return a value or a vector
#' @export
#'
#' @examples
#' l=1
#' r=9
#' s=5
#' LRscore(l,r,s)
LRscore <- function(l,r,s){
  L=l^(1/2)
  R=r^(1/2)
  S=s
  sc=L*R/(S+L*R)
  return(sc)
}