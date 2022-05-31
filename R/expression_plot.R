#' @title Expression Plot
#' @description Displays the level of expression of a gene in each cell on the 2D projected data.
#'
#' @details
#' `name` can be any character that corresponds to a row name
#' of the matrix counts present in the object.
#' @details
#' `colors` must be "default", "rainbow" or "heat" exclusively. "rainbow" and
#' "heat" are the color palettes provided in R.
#'
#' @param obj an SCSRDataModel object
#' @param name the identifier of the gene of interest
#' @param colors "default" returns the default colorpanel, also accepts "rainbow" or "heat"
#' @param most.variables a logical
#'
#' @return The function returns a R plot.
#' @export
#' @importFrom gplots colorpanel
#'
#' @examples
#' print("dataPrepare")
#' data <- matrix(runif(1000,0,1),nrow=50,ncol=20)
#' rownames(data) <- paste("gene",seq_len(50))
#' obj <- dataPrepare(data)
#' print("cell Clustering")
#' obj <- cellClustering(obj)
#' expression_plot(obj, "gene 20")

expression_plot <- function(obj,name,colors=c("default","rainbow","heat"),most.variables=TRUE){

  if (!is(obj, "SCSRDataModel")){
        stop("obj must be a SCSRDataModel object")
    }

  data <- obj@ncounts$matrix
  tsne <- obj@cluster$tsne

  if (!is.null(obj@ncounts$matrix.mv)&most.variables){
        cat("Matrix of most variable genes used. To use the whole matrix set most.variables 
            parameter to FALSE.\n")
        data <- obj@ncounts$matrix.mv
    }

  options(warn=-1)
  colors <- match.arg(colors)
  opar <- par()
  if (is.element(name,rownames(data))==TRUE){
    a <- as.numeric(data[name,])
    a <- a*100/max(a)
    a <- log(5*a+1)
    if (sum(is.na(a))==0){
      if (sum(as.numeric(data[name,]))==0){
        print("No expression")
      }
      if (sum(as.numeric(data[name,]))!=0){
        a <- a+1
        if (colors=="default"){
          cr <- colorpanel(max(a),"lightblue","gray90","#CC0033")
        }
        if (colors=="rainbow"){
          cr <- rainbow(max(a))
        }
        if (colors=="heat"){
          cr <- heat.colors(max(a))
        }
        par(mar=c(5.1, 4.1, 4.1, 8))
        plot(x=tsne[,1],y=tsne[,2],type="n",main=list(name,cex=3, col="black", font=3),xlab="t-SNE1",ylab="t-SNE2")
        bx <- par("usr")
        abline(h=0)
        abline(v=0)
        if (sum(a==1)!=0){
          if (sum(a==1)==1){
            symbols(x=tsne[a==1,1],y=tsne[a==1,2],circles=rep(1,1),inches=0.05,bg=cr[a[a==1]],fg="gray20",add=TRUE)
          } else {
            symbols(x=tsne[a==1,1],y=tsne[a==1,2],circles=rep(1,nrow(tsne[a==1,])),inches=0.05,bg=cr[a[a==1]],fg="gray20",add=TRUE)
          }
        }
        symbols(x=tsne[a!=1,1],y=tsne[a!=1,2],circles=rep(1,nrow(tsne[a!=1,])),inches=0.05,bg=cr[a[a!=1]],fg="gray20",add=TRUE)
        n <- length(cr)
        Y <- bx[3]
        mY <- (bx[4] + bx[3])/2
        dX <- bx[2] - bx[1]
        dY <- bx[4] - bx[3]
        dy <- dY /length(cr)
        dx <- 0.1*dX
        x0 <- bx[2]+dx*0.8
        for (i in seq_len(n)){
          polygon(c(x0,x0+dx,x0+dx,x0), c(Y+(i-1)*dy,Y+(i-1)*dy,Y+i*dy,Y+i*dy), col=cr[i], border=cr[i],xpd=NA)
        }
        mtext("Expression level", side=4,at=mY,xpd=NA,line=+5,cex=1.6)
        mtext("-", side=4,at=Y,xpd=NA,line=+5,cex=1.6,adj=1)
        mtext("+", side=4,at=bx[4],xpd=NA,line=+5,cex=1.6)
        par(opar)
      } else {print("NA not supported")}
    }
  } else {
    print("Change name")
  }
}
