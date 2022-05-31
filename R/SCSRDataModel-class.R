library(methods)

#' SingleCellSignalR Data Model Object
#'
#' An S4 class to represent expression data used for inferring
#' ligand-receptor interactions.
#'
#' @slot ncounts   List containing read count matrix (row names must be set to HUGO 
#' official gene symbols), most variables count matrix if only most variable genes 
#' are kept, vector with genes and parameters used for cleaning.
#' @slot cell.cluster  List containing tsne coordinates, cluster ids and clustering
#' method.
#' @slot cluster  List containing cluster names, number of cells per cluster and
#' markers used to identify each cluster.
#' @slot gene.cluster  List containing list of differentially expressed genes for
#' each cluster and the parameters used for dge.
#' @slot initial.organism         Organism
#'
#' @export
#' @examples
#' new("SCSRDataModel", ncounts = list(matrix(1.5, nrow=2, ncol=2,
#'                                    dimnames=list(c("A","B"), c("C","D"))),
#'                     genes = c("A","B"),
#'                     param = list(normalization = TRUE, outliers = c(0,0))),
#'                     initial.organism = "hsapiens")
#'
setClass("SCSRDataModel",
         slots=c(initial.organism = "character",
                 ncounts = "list",
                 cell.classification = "list",
                 cluster = "list",
                 dge.cluster = "list"
                ),
         prototype = list(
             initial.organism = "hsapiens",
             ncounts = list(matrix = matrix(1, nrow = 2, ncol = 3, dimnames = list(c("A","B"), c("C","D","E"))),
                            matrix.mv = NULL,
                            param = list(outliers = c(0,0), normalization = TRUE, most.variables = 0),
                            genes = c("A","B"), genes.mv = NULL),
             cluster = list(tsne = NULL,
                                id = NULL,
                                method = NULL, 
                                names = NULL),
             cell.classification = list(names = NULL,
                                    id = NULL,
                                    nb.cells = NULL,
                                    markers = NULL),
             dge.cluster = list(genes = NULL,
                                  param = NULL,
                                  markers = NULL)
            
         ))

setValidity("SCSRDataModel",
    function(object) {
        if(!is.character(object@initial.organism))
            return("Initial organism is not a character")
        if(!is.matrix(object@ncounts$matrix))
            return("Specified counts are not a matrix")
        if(!is.numeric(object@ncounts$matrix))
            return("Specified counts are not numeric")
        if (is.null(row.names(object@ncounts$matrix)))
            return("Specified counts have no row names set")
        if (length(object@ncounts$genes) != nrow(object@ncounts$matrix))
            return("Number of genes must match the number of rows in count matrix")

        TRUE
    }
)

setMethod("show", "SCSRDataModel",
    function(object) {
        cat("Organism: ", object@initial.organism, "\n", sep="")
        cat("Expression data:\n")
        if (ncol(object@ncounts$matrix) > 10)
            print(head(object@ncounts$matrix[,1:10]))
        if (!is.null(object@cluster$id)){
            cat("Clusters distribution:\n")
            print(object@cluster$names)
            print(table(object@cluster$id))
        }
        if (!is.null(object@cell.classification$id)){
            cat("Cell classification:\n")
            print(object@cell.classification$nb.cells)
        }
    }
)


# Accessors & setters ========================================================

if (!isGeneric("initialOrganism")) {
    if (is.function("initialOrganism"))
        fun <- initialOrganism
    else
        fun <- function(x) standardGeneric("initialOrganism")
    setGeneric("initialOrganism", fun)
}
#' organism accessor
#' @export
setMethod("initialOrganism", "SCSRDataModel", function(x) x@initial.organism)

if (!isGeneric("ncounts")) {
    if (is.function("ncounts"))
        fun <- ncounts
    else
        fun <- function(x) standardGeneric("ncounts")
    setGeneric("ncounts", fun)
}
#' Counts accessor
#' @export
setMethod("ncounts", "SCSRDataModel", function(x) x@ncounts)

if (!isGeneric("cluster")) {
    if (is.function("cluster"))
        fun <- cluster
    else
        fun <- function(x) standardGeneric("cluster")
    setGeneric("cluster", fun)
}
#' clusters accessor
#' @export
setMethod("cluster", "SCSRDataModel", function(x) x@cluster)

if (!isGeneric("cluster<-")) {
    if (is.function("cluster<-"))
        fun <- `cluster<-`
    else
        fun <- function(x, value) standardGeneric("cluster<-")
    setGeneric("cluster<-", fun)
}
#' clusters setter (internal use only)
setMethod("cluster<-", "SCSRDataModel", function(x,value){
    x@cluster <- value
    methods::validObject(x)
    x
})

if (!isGeneric("cell.classification")) {
    if (is.function("cell.classification"))
        fun <- cell.classification
    else
        fun <- function(x) standardGeneric("cell.classification")
    setGeneric("cell.classification", fun)
}
#' Model cell.classification accessor
#' @export
setMethod("cell.classification", "SCSRDataModel", function(x) x@cell.classification)

if (!isGeneric("cell.classification<-")) {
    if (is.function("cell.classification<-"))
        fun <- `cell.classification<-`
    else
        fun <- function(x, value) standardGeneric("cell.classification<-")
    setGeneric("cell.classification<-", fun)
}
#' cell.cluster setter (internal use only)
setMethod("cell.classification<-", "SCSRDataModel", function(x,value){
    x@cell.classification <- value
    methods::validObject(x)
    x
})

if (!isGeneric("dge.cluster")) {
    if (is.function("dge.cluster"))
        fun <- dge.cluster
    else
        fun <- function(x) standardGeneric("dge.cluster")
    setGeneric("dge.cluster", fun)
}
#' Model cell.cluster accessor
#' @export
setMethod("dge.cluster", "SCSRDataModel", function(x) x@dge.cluster)

if (!isGeneric("dge.cluster<-")) {
    if (is.function("dge.cluster<-"))
        fun <- `dge.cluster<-`
    else
        fun <- function(x, value) standardGeneric("dge.cluster<-")
    setGeneric("dge.cluster<-", fun)
}
#' cell.cluster setter (internal use only)
setMethod("dge.cluster<-", "SCSRDataModel", function(x,value){
    x@dge.cluster <- value
    methods::validObject(x)
    x
})

# Cell clustering ===========================================
if (!isGeneric("cellClustering")) {
    if (is.function("cellClustering"))
        fun <- cellClustering
    else
        fun <- function(obj, ...) standardGeneric("cellClustering")
    setGeneric("cellClustering", fun)
}
#' Clustering
#'
#' Identifies the cell clusters, i.e. the cell subpopulations.
#'
#' @param obj an object of class SCSRDataModel
#' @param n.cluster a number, an estimation of the ideal number of clusters is computed if equal to 0
#' @param n a number, the maximum to consider for an automatic determination of the ideal number of clusters
#' @param method "kmeans" or "simlr"
#' @param plot a logical
#' @param most.variables a logical
#' @param pdf a logical
#' @param write a logical
#'
#' @details If the user knows the number of clusters present in her data set,
#' then `n.cluster` can be set and the estimation of the number of clusters is
#' skipped. `n` is the maximum number of clusters that the automatic estimation
#' of the number of clusters will consider. It is ignored if `n.cluster` is
#' provided. `method` must be "simlr" or "kmeans" exclusively. If set to
#' "simlr", then the function uses the **SIMLR()** function (**SIMLR** package)
#' to perform clustering. If set to "kmeans" the function will perform a
#' dimensionality reduction by principal component analysis (PCA) followed by
#' K-means clustering and 2-dimensional projection by t-distributed stochastic
#' neighbor embedding (t-SNE). Regardless of the value of `method` ("simlr" or
#' "kmeans"), in case `n.cluster` is not provided, then the function relies on
#' the **SIMLR_Estimate_Number_of_Clusters()** function to determine the number
#' of clusters, between 2 and `n`. If `plot` is TRUE, then the function displays
#' the t-SNE map with each cell colored according to the cluster it belongs to.
#' If `method` argument is "simlr", then it further displays a heatmap of the
#' similarity matrix calculated by the **SIMLR()** function. If `pdf` is TRUE,
#' then the function exports the t-SNE plot in a pdf file in the *images*
#' folder. The file is named "t-SNE_map-X.pdf", where X is the `method`
#' argument. If `write` is TRUE, then the function writes two text files in the
#' *data* folder. The first one is called "cluster-Y-X.txt", containing the
#' cluster vector assigning each cell of `data` to a cluster. The second one is
#' called "tsne-Y-X.txt", containing the coordinates of each cell in the 2D
#' t-SNE projection. "X" is the `method` argument anf "Y" is the retained number
#' of clusters. If `most.variables` is TRUE, then the function uses the most 
#' variable genes matrix counts if it exists in the object.
#'
#' @return A SCSRDataModel with cluster definition
#'
#' @export
#'
#' @examples
#' print("dataPrepare")
#' data <- matrix(runif(1000,0,1),nrow=50,ncol=20)
#' rownames(data) <- paste("gene",seq_len(50))
#' obj <- dataPrepare(data)
#' print("cell Clustering")
#' obj <- cellClustering(obj)

setMethod("cellClustering", "SCSRDataModel", function(obj, n.cluster = 0,n = 10, most.variables = TRUE,
                                    method = c("simlr","kmeans"), plot = TRUE, pdf = TRUE,write = TRUE) {
  

    data = obj@ncounts$matrix

    if (!is.null(obj@ncounts$matrix.mv)&most.variables) {
        data = obj@ncounts$matrix.mv
        cat("Matrix of most variable genes used. To use the whole matrix set most.variables 
        parameter to FALSE.\n")
    }

    if (!is.null(obj@cluster$id)){
        cat("Cluster definition will be overwritten.\n")
    }

    if (dir.exists("images")==FALSE & pdf==TRUE){
        dir.create("images")
    }
    if (dir.exists("data")==FALSE & write==TRUE){
      dir.create("data")
    }
    if (n.cluster==0){
        cat("Estimating the number of clusters",fill=TRUE)
        c = SIMLR_Estimate_Number_of_Clusters(data, NUMC=2:n)
        a=data.matrix(as.numeric((c$K1+c$K2)))
        rownames(a)=c(2:n)
        n.cluster = as.numeric(rownames(subset(a, a==min(a))))
        cat(paste("Estimated number of clusters =",n.cluster),fill=TRUE)
    }
    method = match.arg(method)

    final=list()
    if (method=="simlr"){
        s = SIMLR(data,c=n.cluster,no.dim = 2)
        Y = s$y
        cluster = Y$cluster
        final[[1]] = cluster
        final[[2]] = s$F
        final[[3]] = Y$size
        final[[4]] = log(s$S*10^6+1)
        names(final) = c("cluster","t-SNE","numbers","similarity")
    }
    if (method=="kmeans"){
        data = data[rowSums(data)>0,]
        pca = prcomp(t(data),center=TRUE,scale.=TRUE)
        v = (pca$sdev^2)/(sum(pca$sdev^2))
        n = min(which(diff(v)>-10^-4))
        if (n==0 | identical(n, integer(0))){
            n=round(ncol(pca$x)/2)
        }
        pca = pca$x[,seq_len(n)]
        tsne = Rtsne(pca)
        km = kmeans(pca,n.cluster)
        cluster = km$cluster
        final[[1]] = cluster
        final[[2]] = tsne$Y
        final[[3]] = km$size
        names(final) = c("cluster","t-SNE","numbers")
    }

    cr=rainbow(max(cluster))
    if (plot==TRUE){
        if (method=="simlr"){
        pheatmap(log(s$S*10^6+1))
        }
        plot(x=final[[2]][,1],y=final[[2]][,2],type='n',main="t-SNE Map",xlab="t-SNE1",ylab="t-SNE2",xlim=c(min(final[[2]][,1])*1.5,max(final[[2]][,1])*1.1))
        abline(h=0)
        abline(v=0)
        symbols(x=final[[2]][,1],y=final[[2]][,2],circles=rep(1,nrow(final[[2]])),inches=0.04,bg=cr[cluster],add=TRUE)
        legend("topleft",legend = paste("cluster",seq_len(max(cluster))),fill = cr,cex = 0.7)
    }
    if (pdf==TRUE){
        pdf(paste("./images/t-SNE_map-",method,".pdf",sep=""))
        plot(x=final[[2]][,1],y=final[[2]][,2],type='n',main="t-SNE Map",xlab="t-SNE1",ylab="t-SNE2",xlim=c(min(final[[2]][,1])*1.5,max(final[[2]][,1])*1.1))
        abline(h=0)
        abline(v=0)
        symbols(x=final[[2]][,1],y=final[[2]][,2],circles=rep(1,nrow(final[[2]])),inches=0.04,bg=cr[cluster],add=TRUE)
        legend("topleft",legend = paste("cluster",seq_len(max(cluster))),fill = cr,cex = 0.7)
        dev.off()
    }
    if (write==TRUE){
        fwrite(data.frame(final[[1]]),paste("./data/cluster-",n.cluster,"-",method,".txt",sep=""),sep="\t")
        fwrite(data.frame(final[[2]]),paste("./data/tsne-",n.cluster,"-",method,".txt",sep=""),sep="\t")
    }
    cat(paste(n.cluster,"clusters detected"),fill=TRUE)
    for (i in seq_len(n.cluster)){
        cat(paste("cluster",i,"->",final[[3]][i],"cells"),fill=TRUE)
    }

    obj@cluster$method <- method
    obj@cluster$tsne <- final[["t-SNE"]]
    obj@cluster$id <- final$cluster
    obj
})

# Cell classification ===========================================
if (!isGeneric("cellClassifying")){
    if (is.function("cellClassifying"))
        fun <- cellClassifying
    else 
        fun <- function(obj, ...) standardGeneric("cellClassifying")
        setGeneric("cellClassifying",fun)
}

#' Cell classifier
#'
#' Classifies each cellga using cell type specific markers.
#'
#' @param obj an object of type SCSRDataModel
#' @param markers a data frame of cell type signature genes
#' @param plot.details a logical (if TRUE, then plots the number of cells
#' attributed to one cell type, see below)
#' @param most.variables a logical
#' @param write a logical
#' @param verbose a logical
#'
#' @details The ` markers` argument must be a table with cell type gene signatures, one
#' cell type in each column. The column names are the names of the cell types.
#' @details The *markers.default* table provides an example of this format.
#' @details If ` tsne` is not provided in the object, then the function will just not display 
#' the cells on the t-SNE. Although t-SNE maps are widely used to display cells on a 2D
#' projection, the user can provide any table with two columns and a number of
#' rows equal to the number of columns of `data` (e.g. the two first components
#' of a PCA).
#' @details If ` plot.details` is TRUE, then the function plots the number of cells
#' attributed to a single cell type as a function of the threshold applied to
#' the normalized gene signature average.
#' @details If `most.variables` is TRUE, then the function uses the most variable genes
#' matrix counts if it exists in the object.
#' @details If ` write` is TRUE, then the function writes four different text files. (1)
#' The "raw classification matrix" provides the normalized average gene
#' signature for each cell type in each individual cell, a number between 0 and
#' 1. This matrix has one row per cell type and one column per cell, and the
#' sum per column is 1. Row names are the cell type names (column names of the
#' markers table) and the column names are the individual cell identifiers
#' (column names of `data`). (2) The "thresholded classification matrix", which
#' is obtained by eliminating all the values of the "raw classification matrix"
#' that are below a threshold a\*. In practice, a\* is automatically determined
#' by the function to maximize the number of cells that are assigned to a single
#' cell type and all the cells (columns) assigned to 0 or >1 cell types are
#' discarded. The number of cells assigned to a single type depending on a\*
#' can be plotted by using the parameter `plot.details=TRUE`. (3) A cluster
#' vector assigning each cell to a cell type. Note that a supplementary,
#' virtual cluster is created to collect all the cells assigned to 0 or >1
#' types. This virtual cluster is named "undefined". (4) A table associating
#' each cell type to a cluster number in the cluster vector.
#'
#' @return A SCSRDataModel with cluster definition
#'
#' @export
#'
#' @examples
#' print("dataPrepare")
#' data <- matrix(runif(1000,0,1),nrow=50,ncol=20)
#' rownames(data) <- paste("gene",seq_len(50))
#' obj <- dataPrepare(data)
#' print("cell Clustering")
#' obj <- cellClustering(obj)
#'
#' markers <- matrix(paste("gene",seq_len(10)),ncol=5,nrow=2)
#' colnames(markers) <- paste("type",seq_len(5))
#' obj <- cellClassifying(obj, markers = markers)
setMethod("cellClassifying", "SCSRDataModel", function(obj, markers = markers_default,
                                                    plot.details = FALSE, most.variables = TRUE,
                                                    write = TRUE, verbose = TRUE) {
  tsne <- obj@cluster$tsne
  data <- obj@ncounts$matrix
  genes <- obj@ncounts$genes

  if (!is.null(obj@ncounts$matrix.mv)&most.variables){
    cat("Matrix of most variable genes used. To use the whole matrix set most.variables 
        parameter to FALSE.\n")
    data <- obj@ncounts$matrix.mv
    genes <- obj@ncounts$genes.mv
  }

  if (dir.exists("cell-classification")==FALSE & write==TRUE){
      dir.create("cell-classification")
  }
  rownames(data) <- genes
  n.types <- ncol(markers)

  tmp <- data[as.character(unlist(markers))[as.character(unlist(markers)) %in%
                                             genes],]

  if (is.null(dim(tmp))==TRUE){
    stop("Not enough markers genes to pursue the cell classification")
  }
  final <- matrix(0,ncol=ncol(tmp),nrow=(ncol(markers)))

  for (i in seq_len(n.types)){
    m.genes <- markers[,i][markers[,i] %in% genes]
    if (!is.null(dim(tmp[m.genes,]))) final[i,] <- colSums(tmp[m.genes,])/length(m.genes)
    else final[i,] <- NA
  }
  final[is.na(final)] <- 0
  final[,colSums(final)!=0] <- apply(final[,colSums(final)!=0],2,function(x)
    x/sum(x))
  # final <- round(final*1000)/10
  rownames(final) <- colnames(markers)
  colnames(final) <- colnames(data)
  l <- matrix(0,101,2)
  q <- 0
  for (n in seq(0.01,1,0.01)){
    q <- q+1
    f <- final
    f[f<n] <- 0
    l[q+1,1] <- n
    l[q+1,2] <- sum(apply(f,2,function(x) sum(x==0)==(n.types-1)))
  }
  seuil <- max(l[which(l[,2]==max(l[,2])),1])
  m <- matrix(0,(n.types+1),2)
  for (n in 0:n.types){
    f <- final
    f[f<seuil] <- 0
    f <- matrix(f[,apply(f,2,function(x) sum(x==0)==n)],nrow=n.types)
    m[n+1,1] <- n
    m[n+1,2] <- sum(apply(f,2,function(x) sum(x==0)==n))
  }
  m <- matrix(m[m[,2]!=0,],ncol=2)

  res <- list()
  final.s <- final
  final.s[final.s<seuil] <- 0
  res[[1]] <- final.s[,apply(final.s,2,function(x) sum(x==0)==(n.types-1))]
  res[[3]] <- final
  final[!apply(final,2,function(x) x==max(x))] <- 0
  res[[2]] <- final
  h <- pheatmap::pheatmap(res[[1]],cluster_rows=TRUE,cluster_cols=TRUE,
                         show_colnames=FALSE)

  nn <- seq_len(length(rownames(res[[1]])[rowSums(res[[1]])!=0]))
  names(nn) <- rownames(res[[1]])[rowSums(res[[1]])!=0]
  cluster <- unlist(lapply(apply(final.s,2,function(x) names(x[x>0])),
                          function(x) if(length(x)==1){
                            return(nn[as.character(x)])} else {return(0)}))
  cluster[cluster==0]=max(cluster)+1

  if (sum(m[,2]!=0)==1){
    n <- names(nn)
  } else {
    n <- c(names(nn), "Undefined cells")
  }
  cr <- c(rainbow(max(cluster)-1),"gray")
  res[[4]] <- cluster
  # res[[5]] <- h
  res[[5]] <- c(rownames(res[[1]][rowSums(res[[1]])>0,]),"undefined")
  names(res) <- c("tresh_mat","max_mat","raw_mat","cluster","c.names")
  d <- data.frame(Cell_type=c(rownames(res[[1]]),"undefined"),
               Number_of_cells=c(apply(res[[1]],1,function(x) sum(x!=0)),
                                 ncol(res[[3]])-ncol(res[[1]])))
  d <- d[d$Number_of_cells!=0,]
  rownames(d) <- paste("cluster",seq_len(nrow(d)))


  if (verbose==TRUE){
    cat(paste("A threshold of",round(seuil*1000)/10,"% maximizes the number of
              cells assigned to one cell type"),fill=TRUE)
    for (i in (nrow(m)):1){
      cat(paste(m[i,2],"cell(s) or",round(m[i,2]*1000/ncol(data))/10,
                "% identified to",(n.types)-m[i,1],"cell type(s)"),sep=" ",
          fill=TRUE)
    }
    cat(" ",fill=TRUE)
    print(d)
  }

  if (is.null(tsne)==FALSE){
    plot(x=tsne[,1],y=tsne[,2],type='n',main="t-SNE Map",xlab="t-SNE1",
         ylab="t-SNE2",xlim=c(min(tsne[,1])*2,max(tsne[,1])*1.05))
    abline(h=0)
    abline(v=0)
    symbols(x=tsne[,1],y=tsne[,2],circles=rep(1,nrow(tsne)),inches=0.04,
            bg=cr[cluster],add=TRUE)
    legend("topleft",legend=n,fill=cr,cex=0.75)
  }
  if (plot.details==TRUE){
    plot(l,xlab="Threshold (%)",ylab="Number of cells",
         main= "Attribution of one cell to one cell type",type='l')
    abline(v=l[l[,1]==seuil,1], col="red", lty=2)
    symbols(l[l[,1]==seuil,1],l[l[,1]==seuil,2],circles=1,inches=0.05,
            bg="red",fg="red",add=TRUE)
    text(x=l[l[,1]==seuil,1],y=l[l[,1]==seuil,2]-0.07*(max(l)),
         labels=paste(seuil,"%"),offset=20)
  }
  if (write==TRUE){
    fwrite(data.frame(cbind(rownames(res[[1]]),res[[1]])),
           "cell-classification/threshold_class_matrix.txt",sep="\t")
    fwrite(data.frame(cbind(rownames(res[[3]]),res[[3]])),
           "cell-classification/raw_class_matrix.txt",sep="\t")
    nn <- cbind(names(nn),nn)
    colnames(nn) <- c("cluster.name","cluster.number")
    fwrite(data.frame(nn),"cell-classification/cluster2name.txt",sep="\t")
    fwrite(data.frame(cluster),"cell-classification/cluster.txt",sep="\t")
    c.names <- res[[5]]
  }

  obj@cell.classification$id <- res$cluster
  obj@cell.classification$names <- res$c.names
  names(obj@cell.classification$names) = seq(1:length(res$c.names))
  obj@cell.classification$nb.cells <- as.numeric(table(obj@cell.classification$id))
  names(obj@cell.classification$nb.cells) <- obj@cell.classification$names
  obj@cell.classification$markers <- lapply(as.list(markers[,res$c.names[res$c.names!="undefined"]]),
    function(z){ z[!is.na(z) & z != ""]})
  obj
})

# Cluster analysis ===========================================
if (!isGeneric("clusterAnalysis")) {
    if (is.function("clusterAnalysis"))
        fun <- clusterAnalysis
    else
        fun <- function(obj, ...) standardGeneric("clusterAnalysis")
    setGeneric("clusterAnalysis", fun)
}

#' Cluster analysis
#'
#' Analysis of the differentially expressed genes in the clusters
#' and their composition by a marker based approach.
#'
#' @param obj an object of type SCSRDataModel
#' @param dif.exp a logical (if TRUE, then computes the diferential gene
#' expression between the clusters using **edgeR**)
#' @param s.pval a value, a fixed p-value threshold
#' @param markers a table of cell type signature genes
#' @param most.variables a logical
#' @param write a logical
#' @param verbose a logical
#'
#' @details If `dif.exp` is TRUE, then the function uses **edgeR** functions
#' **glmFit()** and **glmRT()** to find differentially expressed genes between
#' one cluster and all the other columns of `data`.
#' @details
#' If `dif.exp` is FALSE, then
#' the function skips the differential gene analysis.
#' @details
#' If `c.names` is not set in the object,
#' the clusters will be named from 1 to the maximum number of
#' clusters (cluster 1, cluster 2, ...). The user can exploit the `c.names`
#' vector in the list returned by the **cell_classifier()** function for this
#' purpose. The user can also provide her own cluster names.
#' @details
#' `s.pval`  is the
#' adjusted (Benjamini-Hochberg) p-value threshold imposed to gene differential
#' expression.
#' @details
#' If `markers` is set, it must be a table with gene signatures for
#' one cell type in each column. The column names are the names of the cell
#' types. If no markers table is set, the function will use the one present in the object.
#' @details
#' If `markers` is not provided, then the function skips the cluster
#' cell type calling step.
#' @details
#' If `write` and `dif.exp` are both TRUE, then the
#' function writes a text file named "table_dge_X.txt", where X is the cluster
#' name, that contains the list of differentially expressed genes.
#' @details 
#' If `most.variables` is TRUE, then the function uses the most variable genes
#' matrix counts if it exists in the object.
#' @details
#' If `write` is TRUE and `markers` is provided, then the function writes in a second text
#' file a table containing probabilities of assignments of each cluster to a
#' cell type for each cell cluster. This cell type calling is performed as for
#' the individual cells without thresholding but based on the cluster average
#' transcriptome.
#'
#' @return A SCSRDataModel with cluster analysis
#'
#' @export
#'
#' @examples
#' print("dataPrepare")
#' data <- matrix(runif(1000,0,1),nrow=50,ncol=20)
#' rownames(data) <- paste("gene",seq_len(50))
#' obj <- dataPrepare(data)
#' print("cell Clustering")
#' obj <- cellClustering(obj)
#' obj <- clusterAnalysis(obj)


setMethod("clusterAnalysis", "SCSRDataModel", function(obj, markers = NULL, dif.exp = TRUE,
                                                    s.pval=10^-2, most.variables = TRUE,
                                                    write = TRUE, verbose = TRUE) {


if (is.null(markers)&!is.null(obj@dge.cluster$markers)){
    markers <- as.data.frame(do.call(cbind, obj@dge.cluster$markers))
    for (c in names(markers)){
        markers[,c][duplicated(markers[,c])] <- NA
    }
}else if (is.null(obj@dge.cluster$markers)){
    obj@dge.cluster$markers <- lapply(as.list(markers),
        function(z){ z[!is.na(z) & z != ""]})
}

cluster <- obj@cluster$id
c.names <- obj@cluster$names
genes <- obj@ncounts$genes
data <- obj@ncounts$matrix

if (!is.null(obj@ncounts$matrix.mv)&most.variables){
    cat("Matrix of most variable genes used. To use the whole matrix set most.variables 
        parameter to FALSE.\n")
    data <- obj@ncounts$matrix.mv
    genes <- obj@ncounts$genes.mv
}

if (dir.exists("cluster-analysis")==FALSE & write==TRUE){
    dir.create("cluster-analysis")
  }
  if (is.null(c.names)==TRUE){
    c.names <- paste("cluster",seq_len(max(cluster)))
  }
  if (min(cluster)!=1){
    cluster <- cluster + 1 - min(cluster)
  }
  if (length(c.names)!=max(cluster) | sum(duplicated(c.names))>0 |
      grepl("/",paste(c.names,collapse =""))){
    stop("The length of c.names must be equal to the number of
        clusters and must contain no duplicates. The cluster names must not
        include special characters")
  }
  rownames(data) <- genes
  n.cluster <- max(cluster)
  z <- c(seq_len(n.cluster))
  class <- list()
  length(class) <- n.cluster
  diff.genes <- list()
  final.2 <- NULL
  types <- NULL
  k <- 0
  v <- NULL
  if (dif.exp==TRUE){
    dge <- DGEList(data,genes=genes)
    dge <- calcNormFactors(dge)
    if (verbose==TRUE){
      cat("edgeR differential gene expression (dge) processing:",fill=TRUE)
    }

    for (i in z){
      if (verbose==TRUE){
        cat(paste("Looking for differentially expressed genes in", c.names[i]),
            fill=TRUE)
      }
      cl <- as.matrix(cluster)
      cl[cluster==i,] <- 1
      cl[cluster!=i,] <- 2
      c <- factor(cl)
      design <- model.matrix(~0+c)
      cm <- makeContrasts(c1-c2,levels=design)
      y <- estimateDisp(dge,design,robust=TRUE)
      fit.y <- glmFit(y,design)
      lrt <- glmLRT(fit.y,contrast=cm)
      sel.r <- topTags(lrt,p.value=s.pval,adjust.method="BH",n=nrow(y$counts))
      if (length(sel.r)==0){
        if (verbose==TRUE){
          cat(paste("No differentially expressed genes in", c.names[i]),
              fill=TRUE)
        }
        v <- c(v,i)
      } else {
        k <- k+1
        resu <- sel.r$table
        resu <- resu[order(resu$logFC,decreasing=TRUE),]
        diff.genes[[k]] <- resu
        if (write==TRUE){
          fwrite(data.frame(resu),paste("./cluster-analysis/table_dge_",
                                        c.names[i],".txt",sep=""),sep="\t")
        }
      }
    }
    if (is.null(v)==TRUE){
      names(diff.genes) <- c.names
    } else {
      names(diff.genes) <- c.names[-v]
    }
  }

  if (is.null(markers)==FALSE){
    final <- matrix(0,nrow=n.cluster,ncol=ncol(markers))
    colnames(final) <- colnames(markers)
    rownames(final) <- c.names
    for (i in z){
      tmp <- NULL
      for (j in seq_len(ncol(markers))){
        m.genes <- markers[,j][markers[,j] %in% genes]
        tmp <- c(tmp,sum(rowSums(data.frame(data[m.genes,
                                                cluster==i]))/sum(cluster==i))/
                  length(m.genes))
        tmp[is.na(tmp)] <- 0
      }
      final[i,] <- as.numeric(tmp)
    }
    final <- final/rowSums(final,na.rm=TRUE)
    final.2 <- final
    m <- apply(final.2,1,function(x) (max(x)-min(x))/2)
    for (i in z){
      final.2[i,][final.2[i,]<m[i]] <- 0
      c.type <- names(which(final.2[i,]!=0))
      if (length(c.type)==1){
        if (verbose==TRUE){
          cat(paste(c.names[i], ": ", c.type, sep=""),fill=TRUE)
        }
        class[[i]] <- c.type
      }
      if (length(c.type)==2){
        if (verbose==TRUE){
          cat(paste(c.names[i], ": ", c.type[1],"/",c.type[2],sep=""),fill=TRUE)
        }
        class[[i]] <- paste(c.type[1],"/",c.type[2],sep="")
      }
      if (length(c.type)>2){
        if (verbose==TRUE){
          cat(paste(c.names[i], ": Unconclusive (more than 2 cell types
                    attributed to the cluster)",sep=""),fill=TRUE)
        }
        class[[i]] <- "Undefined"
      }
    }
    names(class) <- c.names
    pheatmap(t(final.2),cluster_cols=FALSE)
    types <- do.call(rbind,class)
    if (write==TRUE){
      fwrite(data.frame(final),paste("./cluster-analysis/cluster_types.txt",
                                     sep=""),sep="\t")
    }
  }
  res <- list(diff.genes,class,types,final.2,s.pval)
  names(res) <- c("diff.genes","class","types","matrix", "s.pval")

  if (!is.null(markers)){
    if (!is.null(obj@cluster$names)) cat("Cluster names will be overwritten.\n")
    obj@cluster$names <- unlist(class)
    names(obj@cluster$names) <- unlist(strsplit(names(class),"cluster "))[seq_len(
        length(names(class)))*2]
  }
  obj@dge.cluster$genes <- diff.genes
  names(obj@dge.cluster$genes) <- unlist(strsplit(names(diff.genes),"cluster "))[seq_len(
    length(names(diff.genes)))*2]
  obj@dge.cluster$param <- list(s.pval = s.pval)

  obj
})



if (!isGeneric("cellSignaling")) {
    if (is.function("cellSignaling"))
        fun <- cellSignaling
    else
        fun <- function(obj, ...) standardGeneric("cellSignaling")
    setGeneric("cellSignaling", fun)
}

#' Cell Signaling
#'
#' Computes "autocrine" or "paracrine" interactions between cell
#' clusters.
#'
#' @param obj an object of type SCSRDataModel
#' @param int.type "autocrine" or "paracrine"
#' @param s.score LRscore threshold
#' @param logFC a number, the log fold-change threshold for differentially
#' expressed genes
#' @param tol a tolerance parameter for balancing "autocrine|paracrine"
#' interactions to the "autocrine" or "paracrine" group
#' @param most.variables a logical
#' @param write a logical
#' @param verbose a logical
#'
#' @details `int.type` must be equal to "paracrine" or "autocrine" exclusively.
#' The "paracrine" option looks for ligands expressed in cluster A and their
#' associated receptors according to LR*db* that are expressed in any other
#' cluster but A. These interactions are labelled "paracrine". The interactions
#' that involve a ligand and a receptor, both differentially expressed in their
#' respective cell clusters according to the **edgeR** analysis performed by the
#'  **cluster_analysis()** function, are labelled "specific". The "autocrine"
#' option searches for ligands expressed in cell cluster A and their associated
#' receptors also expressed in A. These interactions are labelled "autocrine".
#' Additionally, it searches for those associated receptors in the other cell
#' clusters (not A) to cover the part of the signaling that is "autocrine" and
#' "paracrine" simultaneously. These interactions are labelled
#' "autocrine/paracrine".
#' @details
#' The `tol` argument allows the user to tolerate a fraction of the cells in
#' cluster A to express the receptors in case `int.type="paracrine"`, that is to
#' call interactions that are dominantly paracrine though not exclusively.
#' Conversely, it allows the user to reject interactions involving receptors
#' that would be expressed by a small fraction of cluster A cells in case
#' `int.type="autocrine"`. By construction theassociation of these two options
#' covers all the possible interactions and increasing the `tol` argument allows
#' the user to move interactions from "autocrine" to "paracrine".
#' @details
#' If the user does not set `c.names`, the clusters will be named from 1 to the
#' maximum number of clusters (cluster 1, cluster 2, ...). The user can exploit
#' the `c.names` vector in the list returned by the **cell_classifier()**
#' function for this purpose. The user can also provide her own cluster names.
#' @details
#' `s.score` is the threshold on the LRscore. The value must lie in the [0;1]
#' interval, default is 0.5 to ensure confident ligand-receptor pair
#' identifications (see our publication). Lower values increase the number of
#' putative interactions while increasing the false positives. Higher values
#' do the opposite.
#' @details
#' `logFC` is a threshold applied to the log fold-change (logFC) computed for
#' each gene during the differential gene expression analysis. Its default value
#' is log~2~(1.5) It further selects the differentially expressed genes (>logFC)
#' after the p-value threshold imposed in the function **cluster_analysis()** below.
#' @details
#' `species` must be equal to "homo sapiens" or "mus musculus", default is
#' "homo sapiens". In the case of mouse data, the function converts mouse genes
#' in human orthologs (according to Ensembl) such that LR*db* can be exploited,
#' and finally output genes are converted back to mouse.
#' @details 
#' If `most.variables` is TRUE, then the function uses the most variable genes
#' matrix counts if it exists in the object.
#' @details
#' If `write` is TRUE, then the function writes a text file that reports the
#' interactions in the *cell-signaling* folder. This file is a 4-column table:
#' ligands, receptors, interaction types ("paracrine", "autocrine",
#' "autocrine/paracrine" and "specific"), and the associated LRscore.
#' @details
#' Remarks:
#' @details- This function can be used with any `data` table associated with
#' corresponding `genes` and `cluster` vectors, meaning that advanced users can
#' perform their own data normalization and cell clustering upfront.
#' @details- In case the function **cluster_analysis()** was not executed, this function
#' would work but "specific" interactions would not be annotated as such.
#'
#' @return A SCSRSignaling with LR interaction
#'
#' @export
#'
#' @examples
#' print("dataPrepare")
#' data <- matrix(runif(1000,0,1),nrow=50,ncol=20)
#' rownames(data) <- paste("gene",seq_len(50))
#' obj <- dataPrepare(data)
#' print("cell Clustering")
#' obj <- cellClustering(obj)
#' print("cell Signaling")
#' obj.int <- cellSignaling(obj,int.type = "paracrine")


setMethod("cellSignaling", "SCSRDataModel", function(obj, 
                            int.type = c("paracrine","autocrine"),s.score = 0.5, 
                            logFC = log2(1.5), tol = 0, most.variables = TRUE, write = TRUE, verbose = TRUE) {

    cluster <- obj@cluster$id
    c.names <- obj@cluster$names
    species <- obj@initial.organism
    data <- obj@ncounts$matrix
    param <- list(s.score = s.score, logFC = logFC, tol = tol)

    if (!is.null(obj@ncounts$matrix.mv)&most.variables){
        cat("Matrix of most variable genes used. To use the whole matrix set most.variables 
            parameter to FALSE.\n")
        data <- obj@ncounts$matrix.mv
    }

    if (dir.exists("cell-signaling")==FALSE & write==TRUE){
        dir.create("cell-signaling")
      }
      if (is.null(c.names)==TRUE){
        c.names <- paste("cluster",seq_len(max(cluster)))
      }
      if (min(cluster)!=1){
        cluster <- cluster + 1 - min(cluster)
      }
      if (length(c.names)!=max(cluster) | sum(duplicated(c.names))>0 |
          grepl("/",paste(c.names,collapse =""))){
        stop("The length of c.names must be equal to the number of obj@cluster
            and must contain no duplicates. The obj@cluster names must not include
            special characters")
      }
      int.type <- match.arg(int.type)

      #rownames(data) <- genes
      z <- seq_len(max(cluster))
      lig <- unique(LRdb$ligand)
      rec <- unique(LRdb$receptor)
      data <- data.frame(data)
      data <- data[rowSums(data)>0,]
      med <- sum(data)/(nrow(data)*ncol(data))

      if (species=='mus musculus'){
        Hs2mm <- mm2Hs[,1]
        mm2Hs <- mm2Hs[,2]
        names(mm2Hs) <- as.character(Hs2mm)
        names(Hs2mm) <- as.character(mm2Hs)
        m.names <- mm2Hs[rownames(data)]
        data <- subset(data,(!is.na(m.names)))
        m.names <- m.names[!is.na(m.names)]
        rownames(data) <- as.character(m.names)
      }


      out <- list()
      ## Autocrine -------------------
      if (int.type=="autocrine"){
        if (verbose==TRUE){
          cat("Autocrine signaling: ",fill=TRUE)
        }
        k=0
        int=NULL
        n.int=NULL
        if (verbose==TRUE){
          cat("Checking for cell/cell signaling:",fill=TRUE)
        }
        for (i in z){
          if (sum(cluster==i)>1){
            tmp <- data[,cluster==i]
            tmp <- tmp[rowSums(tmp)>0,]
            if (sum(is.element(lig, rownames(tmp)))>0){
              lig.tmp <- rownames(tmp)[is.element(rownames(tmp),lig)]
              #lig.tmp <- lig.tmp[!is.element(lig.tmp,gene.list[[i]])]
            } else {lig.tmp=NULL}

            final.tmp <- LRdb[is.element(LRdb$ligand,lig.tmp),seq_len(2)]
            final.tmp <- data.frame(final.tmp,as.character(
              rep("autocrine|paracrine",sum(is.element(LRdb$ligand,lig.tmp)))))
            m.lig <- rowSums(tmp[unique(final.tmp[,1]),])/sum(cluster==i)
            names(m.lig) <- unique(final.tmp[,1])

            if (sum(is.element(rec, rownames(tmp)))>0){
              rec.tmp <- rownames(tmp)[is.element(rownames(tmp[apply(
                tmp,1,function(x) sum(x>0))>tol*ncol(tmp),]),rec)]
            } else {rec.tmp=NULL}

            for (j in z){
              if (sum(cluster==i)>1){
                temp <- data[,cluster==j]
                temp <- temp[rowSums(temp)>0,]
                if (sum(is.element(rec, rownames(temp)))>0){
                  rec.temp <- rownames(temp)[is.element(rownames(temp),rec)]
                  #rec.temp <- rec.temp[!is.element(rec.temp,gene.list[[j]])]
                } else {rec.temp=NULL}
                rec.temp <- rec.temp[is.element(rec.temp,rec.tmp)]
                m.rec <- rowSums(data.frame(temp[rec.temp,]))/sum(cluster==j)
                names(m.rec) <- rec.temp

                final <- final.tmp[is.element(final.tmp$receptor,rec.temp),]
                final <- cbind(final,LRscore(m.lig[final$ligand],m.rec[final$receptor],
                                          med))

                colnames(final) <- c(c.names[i],c.names[j],"interaction type","LRscore")

                if (i==j){
                  final$`interaction type`="autocrine"
                }

                final <- final[final[,4]>s.score,]
                final <- final[order(final[,4],decreasing=TRUE),]

                if (species=="mus musculus"){
                  final[,1] <- Hs2mm[as.character(final[,1])]
                  final[,2] <- Hs2mm[as.character(final[,2])]
                }

                if (nrow(final)>0){
                  k <- k+1
                  out[[k]] <- final
                  if (verbose==TRUE){
                    cat(paste(nrow(final),"interactions from",c.names[i],
                              "to",c.names[j]),fill=TRUE)
                  }
                  int <- c(int,paste(i,"-",j,sep=""))
                  n.int <- c(n.int,paste(c.names[i],"-",c.names[j],sep=""))
                  gr <- graph_from_data_frame(final,directed=FALSE)
                  if (write==TRUE){
                    fwrite(data.frame(final),paste("./cell-signaling/LR_interactions_",
                                                   c.names[i],"-",c.names[j],"-",
                                                   int.type,".txt",sep=""),sep="\t")
                  }
                }
              }
            }
          }
        }
        if (k!=0){
          names(out) <- n.int
        }
      }


      ## Paracrine -------------------
      if (int.type=="paracrine"){
        gene.list <- vector("list",max(cluster))
        for (i in z){
          if (is.null(obj@dge.cluster$genes)&
                !file.exists(paste("./cluster-analysis/table_dge_",c.names[i],
                                ".txt",sep=""))) {
            gene.list[[i]] <- "none"
            cat(paste("No such file as table_dge_",c.names[i],
                      ".txt in the cluster-analysis folder", sep=""),fill=TRUE)
          }else{
            if (file.exists(paste("./cluster-analysis/table_dge_",c.names[i],
                                ".txt",sep=""))) {
                resu <- fread(paste("./cluster-analysis/table_dge_",c.names[i],
                               ".txt",sep=""),data.table=FALSE)
            }else{
                resu <- obj@dge.cluster$genes[[i]]
            }
            gene.list[[i]] <- resu$genes[resu$logFC>logFC]
            if (species == "mus musculus"){
              gene.list[[i]] <- mm2Hs[gene.list[[i]]]
            }
            gene.list[[i]] <- gene.list[[i]][!is.na(gene.list[[i]])]
          }
        }

        if (verbose==TRUE){
          cat("Paracrine signaling: ",fill=TRUE)
        }
        k=0
        int=NULL
        n.int=NULL
        if (verbose==TRUE){
          cat("Checking for signaling between cell types",fill=TRUE)
        }
        for (i in z){
          if (sum(cluster==i)>1){
            tmp <- data[,cluster==i]
            tmp <- tmp[rowSums(tmp)>0,]
            if (sum(is.element(lig, rownames(tmp)))>0){
              lig.tmp <- rownames(tmp)[is.element(rownames(tmp),lig)]
              #lig.tmp <- lig.tmp[!is.element(lig.tmp,gene.list[[i]])]
            } else {lig.tmp=NULL}

            final.tmp <- LRdb[is.element(LRdb$ligand,lig.tmp),seq_len(2)]
            final.tmp <- data.frame(final.tmp,as.character(
              rep("paracrine",sum(is.element(LRdb$ligand,lig.tmp)))),
              stringsAsFactors=FALSE)
            m.lig <- rowSums(tmp[unique(final.tmp[,1]),])/sum(cluster==i)
            names(m.lig) <- unique(final.tmp[,1])

            if (sum(is.element(rec, rownames(tmp)))>0){
              rec.tmp <- rownames(tmp)[is.element(rownames(tmp[
                apply(tmp,1,function(x) sum(x>0))>tol*ncol(tmp),]),rec)]
            } else {rec.tmp=NULL}

            for (j in z[-i]){
              if (sum(cluster==j)>1){
                temp <- data[,cluster==j]
                temp <- temp[rowSums(temp)>0,]
                if (sum(is.element(rec, rownames(temp)))>0){
                  rec.temp <- rownames(temp)[is.element(rownames(temp),rec)]
                  #rec.temp <- rec.temp[!is.element(rec.temp,gene.list[[j]])]
                } else {rec.temp=NULL}
                rec.temp <- rec.temp[!is.element(rec.temp,rec.tmp)]
                m.rec <- rowSums(data.frame(temp[rec.temp,]))/sum(cluster==j)
                names(m.rec) <- rec.temp

                final <- final.tmp[is.element(final.tmp$receptor,rec.temp),]
                final <- cbind(final,LRscore(m.lig[final$ligand],m.rec[final$receptor],
                                          med))
                exclus <- final$ligand %in% gene.list[[i]] & final$receptor %in%
                  gene.list[[j]]
                if (sum(exclus)!=0){
                  f.exclu <- final[exclus,]
                  final <- final[!(final$ligand %in% gene.list[[i]] & final$receptor
                                  %in% gene.list[[j]]),]
                  f.exclu[,3] <- "specific"
                  final <- rbind(f.exclu,final)
                }

                colnames(final) <- c(c.names[i],c.names[j],"interaction type",
                                    "LRscore")
                final <- final[final[,4]>s.score,]
                final <- final[order(final[,4],decreasing=TRUE),]

                if (species=="mus musculus"){
                  final[,1] <- Hs2mm[as.character(final[,1])]
                  final[,2] <- Hs2mm[as.character(final[,2])]
                }

                if (nrow(final)>0){
                  k=k+1
                  out[[k]] <- final
                  if (verbose==TRUE){
                    cat(paste(nrow(final),"interactions from",c.names[i],"to",
                              c.names[j]),fill=TRUE)
                  }
                  int <- c(int,paste(i,"-",j,sep=""))
                  n.int <- c(n.int,paste(c.names[i],"-",c.names[j],sep=""))
                  gr <- graph_from_data_frame(final,directed=FALSE)
                  if (write==TRUE){
                    fwrite(data.frame(final),paste(
                      "./cell-signaling/LR_interactions_",c.names[i],"-",c.names[j],
                      "-",int.type,".txt",sep=""),sep="\t")
                  }
                } else {
                  if (verbose==TRUE){
                    cat(paste(nrow(final),"No significant interaction found from",c.names[i],"to",
                                                          c.names[j]),fill=TRUE)
                  }
                }
              }
            }
          }
        }
        if (k!=0){
          names(out) <- n.int
        }
    }
    if (length(out)>0){
        out_subpop <- data.frame(Subpop=names(out))
        out_subpop$Number <- unlist(lapply(out, nrow))
        LR <- foreach (inter = 1:length(out), .combine = 'rbind') %do%{
            data.frame(Ligand = out[[inter]][,1], Receptor = out[[inter]][,2])
            }
        new("SCSRInteraction", LRinter = out, LRsubpop = out_subpop, ligands = LR$Ligand,
        receptors = LR$Receptor, param = param)
    }else{
        new("SCSRInteraction", LRinter = NULL, LRsubpop = NULL, ligands = NULL,
        receptors = NULL, param = param)
    }
}) 
