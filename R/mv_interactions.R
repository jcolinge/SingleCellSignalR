#' @title most variable interactions
#' @description Displays a heatmap showing the most variable interactions over
#' all clusters.
#'
#'
#' @param obj an object of type SCSRDataModel
#' @param n an integer the number of most variables interactions
#' @param most.variables a logical
#'
#' @return The function displays a heatmap showing the most variable
#' interactions over all clusters
#' @export
#'
#' @importFrom pheatmap pheatmap
#' @importFrom stats var
#'
#' @examples
#'
#' print("dataPrepare")
#' data <- matrix(runif(1000,0,1),nrow=50,ncol=20)
#' rownames(data) <- paste("gene",seq_len(50))
#' obj <- dataPrepare(data, cluster = c(rep(1,10),rep(2,10)))

#' mv_interactions(obj)

mv_interactions <- function(obj,n=30,most.variables=TRUE){


  if (!is(obj, "SCSRDataModel"))
        stop("obj must be an object of class SCSRDataModel")

  c.names <- obj@cluster$names
  cluster <- obj@cluster$id
  genes <- obj@ncounts$genes
  data <- obj@ncounts$matrix
   if (!is.null(obj@ncounts$matrix.mv)&most.variables){
        cat("Most variable genes used. To use the whole gene set set most.variables 
            parameter to FALSE.\n")
        data <- obj@ncounts$matrix.mv
        genes <- obj@ncounts$genes.mv
    }
  species <- obj@initial.organism

  if (is.null(c.names)==TRUE){
    c.names <- paste("cluster",seq_len(max(cluster)))
  }

  if (species=='mus musculus'){
    Hs2mm <- mm2Hs[,1]
    mm2Hs <- mm2Hs[,2]
    names(mm2Hs) <- Hs2mm
    names(Hs2mm) <- as.character(mm2Hs)
    m.names <- mm2Hs[genes]
    data <- subset(data,(!is.na(m.names)))
    m.names <- m.names[!is.na(m.names)]
    genes <- as.character(m.names)
  }
  rownames(data) <- genes
  l_sc <- LRdb[is.element(LRdb$ligand,rownames(data)),]
  int_sc <- l_sc[is.element(l_sc$receptor,rownames(data)),]
  lr_sc <- matrix(0,nrow=nrow(int_sc),ncol=max(cluster)^2)
  rownames(lr_sc) <- paste(int_sc$ligand,int_sc$receptor,sep=" / ")
  med <- sum(data)/(nrow(data)*ncol(data))
  nam <- NULL
  q <- 0
  for (i in seq_len(max(cluster))){
    for (j in seq_len(max(cluster))){
      q <- q+1
      lr_sc[,q] <- (rowMeans(data[int_sc$ligand,cluster==i])*
                     rowMeans(data[int_sc$receptor,cluster==j]))^0.5/
        (med + (rowMeans(data[int_sc$ligand,cluster==i])*rowMeans(
          data[int_sc$receptor,cluster==j]))^0.5)
      nam=c(nam,paste(c.names[i],c.names[j],sep=" -> "))
    }
  }
  colnames(lr_sc) <- nam
  if (sum(lr_sc)!=0){
    if (nrow(lr_sc)<n){
      n <- nrow(lr_sc)
    }
    lr_sc <- subset(lr_sc,rowSums(lr_sc)!=0)
    v <- apply(lr_sc,1,var)/apply(lr_sc,1,mean)
    lr_sc <- lr_sc[order(v,decreasing=TRUE),]
    lr_sc <- lr_sc[apply(lr_sc,1, max)>0.5,]
    pheatmap::pheatmap(lr_sc[seq_len(n),colSums(lr_sc[seq_len(n),])!=0],
                       cluster_cols=TRUE)
  } else {
    cat("No interactions detected. Make sure the genes vector is composed of
        HUGO official gene names.",fill=TRUE)
  }

}


