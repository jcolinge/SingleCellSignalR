#' @title most variable interactions
#' @description Displays a heatmap showing the most variable interactions over
#' all clusters.
#'
#'
#' @param obj an object of type SCSRDataModel
#' @param n an integer the number of most variables interactions
#' @param most.variables a logical
#' @param addLR a dataframe with 2 columns
#'
#' @details 
#' The `addLR` allows the user to add LR interaction that they want to study.
#' It must be a dataframe with one column "ligand" and one column "receptor".
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

mv_interactions <- function(obj,n=30,most.variables=TRUE,addLR = NULL){



  if (!is(obj, "SCSRDataModel"))
       { stop("obj must be an object of class SCSRDataModel")}

  c.names <- obj@cluster$names
  cluster <- obj@cluster$id
  data <- obj@ncounts$matrix
  species <- obj@initial.organism
  

  if (is.null(c.names)){
    c.names <- paste("cluster",seq_len(max(cluster)))
  }
  
  if (!is.null(addLR)){
      if (!is.data.frame(addLR)|ncol(addLR)<2){
          cat("Please input the ligand-receptors interactions you want to add as
                a two-column dataframe.\n")
        }else {
          if (!any(addLR[,1]%in%genes)){
              if (species != "hsapiens"){
                    ortho <- data.frame(Hsapiens = rownames(object@ncounts$matrix),
                                        species = object@ncounts$initial.orthologs)
                    lig <- data.frame(Hsapiens = addLR[,1], id = seq(1, nrow(addLR), 1))
                    rec <- data.frame(Hsapiens = addLR[,2], id = seq(1, nrow(addLR), 1))
                    lig <- merge(lig, ortho, by.x = "Hsapiens", order = FALSE)
                    lig <- lig[order(lig$id),]
                    rec <- merge(rec, ortho, by.x = "Hsapiens", order = FALSE)
                    rec <- rec[order(rec$id),]
                    addLR[,1] <- lig$species
                    addLR[,2] <- rec$species
                }
          }
          if (species != "hsapiens"){
              ortho <- data.frame(Hsapiens = rownames(object@ncounts$matrix),
                                        species = object@ncounts$initial.orthologs)
              lig <- data.frame(Hsapiens = LRdb$ligand, id = seq(1, nrow(LRdb), 1))
              rec <- data.frame(Hsapiens = LRdb$receptor, id = seq(1, nrow(LRdb), 1))
              lig <- merge(lig, ortho, by.x = "Hsapiens", order = FALSE)
              lig <- lig[order(lig$id),]
              rec <- merge(rec, ortho, by.x = "Hsapiens", order = FALSE)
              rec <- rec[order(rec$id),]
              LRdb$ligand <- lig$species
              LRdb$receptor <- rec$species
          }

          add <- cbind(addLR[,1:2],data.frame(matrix(NA,
                                    ncol=ncol(LRdb)-2),nrow=nrow(addLR)))
          names(add) <- names(LRdb)
          LRdb <- rbind(LRdb, add)
          LRdb <- unique(LRdb)
      }
    }
    if (species!='hsapiens'){
    rownames(data) <- obj@ncounts$initial.orthologs
  } 
   if (!is.null(obj@ncounts$matrix.mv)&most.variables){
        cat("Most variable genes used. To use the whole gene set set most.variables 
            parameter to FALSE.\n")
        data <- obj@ncounts$matrix.mv
        if (species!='hsapiens'){
          rownames(data) <- obj@ncounts$initial.orthologs.mv
        } 
    }
  
  genes <- rownames(data)

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
      if (length(data[int_sc$ligand,cluster==i])>1){
        lr_sc[,q] <- (rowMeans(data[int_sc$ligand,cluster==i])*
          rowMeans(data[int_sc$ligand,cluster==j]))^0.5/
        (med + (rowMeans(data[int_sc$ligand,cluster==i])*
          rowMeans(data[int_sc$ligand,cluster==j]))^0.5)
      }else{
        lr_sc[,q] <- ((data[int_sc$ligand,cluster==i])*
          (data[int_sc$ligand,cluster==j]))^0.5/
        (med + ((data[int_sc$ligand,cluster==i])*
          (data[int_sc$ligand,cluster==j]))^0.5)
      }
    
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
                       cluster_cols=TRUE, cellheight=8, cellwidth = 10)
  } else {
    cat("No interactions detected. Make sure the genes vector is composed of
        HUGO official gene names.",fill=TRUE)
  }

}


