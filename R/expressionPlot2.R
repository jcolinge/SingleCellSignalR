#' @title Expression Plot 2
#' @description Displays the level of expression of two genes in each cell on the 2D projected data.
#'
#' @details This function can be used independantly from any other. It displays
#' the expression level of two genes of interest on a 2D projection.
#' @details
#' `name.1` and `name.2` can be any characters that correspond to a row name of the count
#' matrix present in the object.
#'
#' @param obj an SCSRDataModel object
#' @param name.1 the identifier of the first gene of interest
#' @param name.2 the identifier of the second gene of interest
#' @param most.variables a logical
#'
#' @return The function returns a R plot.
#' @export
#'
#' @examples
#' print("dataPrepare")
#' data <- matrix(runif(1000,0,1),nrow=50,ncol=20)
#' rownames(data) <- paste("gene",seq_len(50))
#' obj <- dataPrepare(data)
#' print("cell Clustering")
#' obj <- cellClustering(obj)
#' expression_plot(obj, "gene 20", "gene 2")
#'
expression_plot_2 <- function(obj,name.1,name.2,most.variables=TRUE){
  if (!is(obj, "SCSRDataModel")){
        stop("obj must be a SCSRDataModel object")
    }

  data <- obj@ncounts$matrix
  tsne <- obj@cell.representation$coordinates
  species <- obj@initial.organism
  if (species!='hsapiens'){
    rownames(data) <- obj@ncounts$initial.orthologs
  } 

  if (!is.null(obj@ncounts$matrix.mv)&most.variables){
        cat("Matrix of most variable genes used. To use the whole matrix set most.variables 
            parameter to FALSE.\n")
        data <- obj@ncounts$matrix.mv
        if (species!='hsapiens'){
          rownames(data) <- obj@ncounts$initial.orthologs.mv
        } 
    }

  if (is.element(name.1,rownames(data))==TRUE & is.element(name.2,rownames(data))==TRUE){
    options(warn=-1)
    a <- as.numeric(data[name.1,])
    a <- a*100/max(a)
    a <- log(5*a+1)
    b <- as.numeric(data[name.2,])
    b <- b*100/max(b)
    b <- log(5*b+1)
    opar <- par()
    if (sum(is.na(a))==0){
      if (sum(as.numeric(data[name.1,]))==0){
        print(paste(name.1,"-> No expression"))
      }
      if (sum(as.numeric(data[name.2,]))==0){
        print(paste(name.2,"-> No expression"))
      }
      if (sum(as.numeric(data[name.1,]))!=0 & sum(as.numeric(data[name.2,]))!=0){
        a <- a+1
        b <- b+1
        cr.1 <- colorRampPalette(c("ivory","gray90","#CC0033"),alpha=TRUE)(max(a))
        cr.2 <- colorRampPalette(c("ivory","gray90","dodgerblue3"),alpha=TRUE)(max(b))
        par(mar=c(5.1, 4.1, 4.1, 8))
        plot(x=tsne[,1],y=tsne[,2],type="n",main=list(paste(name.1,"/",name.2),cex=3, col="black", font=3),xlab="t-SNE1",ylab="t-SNE2")
        bx <- par("usr")
        abline(h=0)
        abline(v=0)
        if (sum(a==1)!=0){
          symbols(x=tsne[a==1,1],y=tsne[a==1,2],circles=rep(1,nrow(tsne[a==1,])),inches=0.05,bg=cr.1[a[a==1]],fg="gray20",add=TRUE)
        }
        if (sum(b==1)!=0){
          symbols(x=tsne[b==1,1],y=tsne[b==1,2],circles=rep(1,nrow(tsne[b==1,])),inches=0.05,bg=cr.2[b[b==1]],fg="gray20",add=TRUE)
        }
        symbols(x=tsne[a!=1,1],y=tsne[a!=1,2],circles=rep(1,nrow(tsne[a!=1,])),inches=0.05,bg=cr.1[a[a!=1]],fg="gray20",add=TRUE)
        symbols(x=tsne[b!=1,1],y=tsne[b!=1,2],circles=rep(1,nrow(tsne[b!=1,])),inches=0.05,bg=cr.2[b[b!=1]],fg="gray20",add=TRUE)
        n.1 <- length(cr.1)
        DY <- (bx[4]-bx[3])/2
        Y <- bx[3]
        dX <- bx[2] - bx[1]
        dY <- (bx[4] - bx[3])/2
        dy <- dY / n.1
        dx <- 0.1*dX
        x0 <- bx[2]+dx*0.8
        for (i in seq_len(n.1)){
          polygon(c(x0,x0+dx,x0+dx,x0), c(Y+(i-1)*dy,Y+(i-1)*dy,Y+i*dy,Y+i*dy), col=cr.1[i], border=cr.1[i],xpd=NA)
        }
        mtext(name.1, side=4,at=Y+(n.1/2)*dy,xpd=NA,line=+5,cex=1.2)
        mtext("-", side=4,at=Y+1,xpd=NA,line=+5,cex=1.6,adj=1)
        mtext("+", side=4,at=Y+(n.1)*dy-1,xpd=NA,line=+5,cex=1.6)

        n.2 <- length(cr.2)
        Y <- (bx[4]+bx[3])/2+0.1*bx[4]
        dy <- dY / n.2
        for (i in seq_len(n.2)){
          polygon(c(x0,x0+dx,x0+dx,x0), c(Y+(i-1)*dy,Y+(i-1)*dy,Y+i*dy,Y+i*dy), col=cr.2[i], border=cr.2[i],xpd=NA)
        }
        mtext(name.2, side=4,at=Y+(n.2/2)*dy,xpd=NA,line=+5,cex=1.2)
        mtext("-", side=4,at=Y+1,xpd=NA,line=+5,cex=1.6,adj=1)
        mtext("+", side=4,at=Y+(n.2)*dy-1,xpd=NA,line=+5,cex=1.6)
        par(opar)
      } else {print("NA not supported")}
    }
  } else {
    print("Change names")
  }
}
