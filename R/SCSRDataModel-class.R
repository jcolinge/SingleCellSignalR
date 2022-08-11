#' SingleCellSignalR Data Model Object
#'
#' An S4 class to represent expression data used for inferring
#' ligand-receptor interactions.
#'
#' @slot initial.organism         Organism
#' @slot ncounts   List containing read count matrix (row names must be   
#' set to HUGO official gene symbols), most variables count matrix if    
#' only most variable genes are kept, vector with genes and parameters   
#' used for cleaning.
#' @slot cluster  List containing cluster names, number of cells per  
#' cluster and markers used to identify each cluster.
#' @slot cell.representation  List containing tsne coordinates, cluster  
#' ids and clustering method.
#' @slot dge.cluster  List containing list of differentially expressed  
#' genes for each cluster and the parameters used for dge.
#'
#' @import methods
#'
#' @export
#' @examples
#' new('SCSRDataModel', ncounts = list(matrix=matrix(1, nrow=2, ncol=2,
#'                                    dimnames=list(c('A','B'), c('C','D'))),
#'                     genes = c('A','B'),
#'                     param = list(normalization = TRUE, outliers = c(0,0))),
#'                     initial.organism = 'hsapiens')
#'
setClass("SCSRDataModel",
    slots = c(initial.organism = "character",
        ncounts = "list",
        cluster = "list",
        cell.representation = "list",
        dge.cluster = "list"),
    prototype = list(initial.organism = "hsapiens",
        ncounts = list(matrix = matrix(1,
            nrow = 2, ncol = 3,
            dimnames = list(c("A",
                "B"), c("C",
                "D", "E"))),
            matrix.mv = NULL,
            initial.orthologs = NULL,
            initial.orthologs.mv = NULL,
            param = list(formating = "none",
                specific = list(),
                most.variables = 0)),
        cluster = list(id = NULL,
            method = NULL,
            names = NULL,
            markers = NULL),
        cell.representation = list(coordinates = NULL,
            method = NULL),
        dge.cluster = list(genes = NULL,
            param = NULL)))

setValidity("SCSRDataModel", function(object) {
    if (!is.character(object@initial.organism))
        return("Initial organism is not a character")
    if (!is.matrix(object@ncounts$matrix))
        return("Specified counts are not a matrix")
    if (!is.numeric(object@ncounts$matrix))
        return("Specified counts are not numeric")
    if (is.null(row.names(object@ncounts$matrix)))
        return("Specified counts have no row names set")

    TRUE
})

setMethod("show", "SCSRDataModel", function(object) {
    cat("Organism: ", object@initial.organism, "\n",
        sep = "")
    cat("Expression data:\n")
    if (ncol(object@ncounts$matrix) > 10)
        print(head(object@ncounts$matrix[, seq_len(10)]))
    if (!is.null(object@cluster$id)) {
        cat("Clusters distribution:\n")
        print(object@cluster$names)
        print(table(object@cluster$id))
    }
})


# Accessors & setters
# ========================================================

if (!isGeneric("initialOrganism")) {
    if (is.function("initialOrganism")) {
        fun <- initialOrganism
    } else fun <- function(x) standardGeneric("initialOrganism")

    setGeneric("initialOrganism", fun)
}
#' Generic accessor for the organism
#' @name initialOrganism
#' @aliases initialOrganism,SCSRDataModel-method
#'
#' @param x A SCSRDataModel object
#'
#' @export
#' @keywords internal
#' 
#' @examples
#' message('initialOrganism')
setMethod("initialOrganism", "SCSRDataModel",
    function(x) x@initial.organism)

if (!isGeneric("ncounts")) {
    if (is.function("ncounts"))
        fun <- ncounts 
    else fun <- function(x) standardGeneric("ncounts")

    setGeneric("ncounts", fun)
}
#' Generic accessor for the ncounts slot
#' @name ncounts
#' @aliases ncounts,SCSRDataModel-method
#'
#' @param x A SCSRDataModel object
#'
#' @export
#' @keywords internal
#' 
#' @examples
#' message('ncounts')
setMethod("ncounts", "SCSRDataModel", function(x) x@ncounts)

if (!isGeneric("cluster")) {
    if (is.function("cluster"))
        fun <- cluster 
    else fun <- function(x) standardGeneric("cluster")
    setGeneric("cluster", fun)
}
#' Generic accessor for the cluster slot
#' @name cluster
#' @aliases cluster,SCSRDataModel-method
#'
#' @param x A SCSRDataModel object
#'
#' @export
#' @keywords internal
#' 
#' @examples
#' message('cluster')
setMethod("cluster", "SCSRDataModel", function(x) x@cluster)

if (!isGeneric("cluster<-")) {
    if (is.function("cluster<-"))
        fun <- `cluster<-` 
    else fun <- function(x, value) standardGeneric("cluster<-")
    setGeneric("cluster<-",
        fun)
}
#' clusters setter (internal use only)
#' @param x A SCSRDataModel object
#' @param value valut to be set for SCSRDataModel
#' @keywords internal
setMethod("cluster<-", "SCSRDataModel", function(x, value) {
    x@cluster <- value
    methods::validObject(x)
    x
})

if (!isGeneric("cell.representation")) {
    if (is.function("cell.representation"))
        fun <- cell.representation 
    else fun <- function(x) standardGeneric("cell.representation")
    setGeneric("cell.representation", fun)
}
#' Generic accessor for the cell.representation slot
#' @name cell.representation
#' @aliases cell.representation,SCSRDataModel-method
#'
#' @param x A SCSRDataModel object
#'
#' @export
#' @keywords internal
#' 
#' @examples
#' message('cell.representation')
setMethod("cell.representation", "SCSRDataModel",
    function(x) x@cell.representation)

if (!isGeneric("cell.representation<-")) {
    if (is.function("cell.representation<-"))
        fun <- `cell.representation<-` 
    else fun <- function(x, value) standardGeneric("cell.representation<-")
    setGeneric("cell.representation<-",
        fun)
}
#' cell.representation setter (internal use only)
#' @param x A SCSRDataModel object
#' @param value valut to be set for SCSRDataModel
#' @keywords internal
setMethod("cell.representation<-", "SCSRDataModel",
    function(x, value) {
        x@cell.representation <- value
        methods::validObject(x)
        x
    })

if (!isGeneric("dge.cluster")) {
    if (is.function("dge.cluster"))
        fun <- dge.cluster 
    else fun <- function(x) standardGeneric("dge.cluster")
    setGeneric("dge.cluster", fun)
}
#' Generic accessor for the dge.cluster slot
#' @name dge.cluster
#' @aliases dge.cluster,SCSRDataModel-method
#'
#' @param x A SCSRDataModel object
#'
#' @export
#' @keywords internal
#' 
#' @examples
#' message('dge.cluster')
setMethod("dge.cluster", "SCSRDataModel",
    function(x) x@dge.cluster)

if (!isGeneric("dge.cluster<-")) {
    if (is.function("dge.cluster<-"))
        fun <- `dge.cluster<-` 
    else fun <- function(x, value) standardGeneric("dge.cluster<-")
    setGeneric("dge.cluster<-",
        fun)
}
#' dge.cluster setter (internal use only)
#' @param x A SCSRDataModel object
#' @param value valut to be set for SCSRDataModel
#' @keywords internal
setMethod("dge.cluster<-", "SCSRDataModel", function(x, value) {
    x@dge.cluster <- value
    methods::validObject(x)
    x
})

# Add cell clustering
# ===========================================
if (!isGeneric("addClustering")) {
    if (is.function("addClustering"))
        fun <- addClustering 
    else fun <- function(obj, ...) standardGeneric("addClustering")
    setGeneric("addClustering", fun)
}
#' addClustering
#'
#' Add the cell clusters, id, names and 2D coordinates for graphs
#' @name addClustering
#' @aliases addClustering,SCSRDataModel-method
#'
#' @param obj an object of class SCSRDataModel
#' @param cluster.id a dataframe with cluster id and names
#' @param coordinates a dataframe with x and y coordinates for the projection
#' @param method a character, giving the projection method
#' @param plot a logical
#'
#' @details The ` obj` argument must be a a SCSRDataModel object, obtained 
#' by applying the ` dataPrepare()` function to the counts matrix (either raw 
#' or already normalized).
#' @details The ` cluster.id` argument must be a dataframe containing 2 
#' columns: cluster  ids and names. It must have as many rows as there are 
#' columns in the matrix counts.
#' @details The ` coordinates` argument must be a dataframe with 2 columns 
#' containing x and y coordinates in the 2D projection. It must have as many 
#' rows as there are columns in the matrix counts.
#' @details The ` method` defines the type of 2D projection (umap, pca, 
#' tsne, ...).
#' @details The ` plot` argument must be a logical to chose whether to plot 
#' the projection
#' @details The ` coordinates` and ` method` arguments are not mandatory but 
#' are needed in order to do cell visualization.
#'
#' @return A SCSRDataModel with cluster definition
#'
#' @export
#'
#' @examples
#' message('addClustering')
#' message('--dataPrepare')
#' data <- matrix(runif(1000,0,1),nrow=50,ncol=20)
#' rownames(data) <- paste('gene',seq_len(50))
#' obj <- dataPrepare(data)
#' message('--add Clustering')
#' obj <- addClustering(obj,cluster.id = sample(1:5,size = 20,replace = TRUE))

setMethod("addClustering", "SCSRDataModel", function(obj, cluster.id, 
    coordinates = NULL, method = NULL, plot = TRUE) {


    if (!is.data.frame(cluster.id)) {
        cluster.id <- as.data.frame(cluster.id)
        names(cluster.id) <- "id"
    }
    if (is.factor(cluster.id[, 1])) {
        cluster.id[, 1] <- as.numeric(cluster.id[, 1])
    }
    if (ncol(obj@ncounts$matrix) != nrow(cluster.id)) {
        message("Cluster.id must contain as many elements as there are cells
            (columns) in the object matrix counts.")
    } else if (!is.null(coordinates) & is.null(method)) {
        message("Please input the projection method corresponding to the 
            coordinates.")
    } else {
        if (!is.null(coordinates)) {
            if (ncol(obj@ncounts$matrix) != nrow(coordinates)) {
                message("Coordinates must contain as many elements as there 
                    are cells (columns) in the object matrix counts.")
            } else if (ncol(coordinates) != 2) {
                message("Coordinates must be a dataframe with 2 columns: x 
                    and y.")
            } else {
                coordinates <- as.data.frame(coordinates)
                names(coordinates) <- c("x", "y")
                obj@cell.representation$coordinates <- coordinates
                obj@cell.representation$method <- method
            }
        } else {
            message("No coordinates given, plotting of cell projection will  
                not be possible.")
        }
        if (ncol(cluster.id) != 2) {
            message("No cluster names inputted.")
            cluster.id$names <- cluster.id[, 1]
        }
        obj@cluster$id <- cluster.id[, 1]
        names <- unique(cluster.id)
        names <- names[order(names[, 1]), ]
        names.vec <- names[, 2]
        names(names.vec) <- names[, 1]

        obj@cluster$names <- names.vec
    }

    if (plot & !is.null(method) & !is.null(coordinates)) {
        cr <- rainbow(max(cluster.id[, 1]))
        plot(x = coordinates$x, y = coordinates$y, type = "n",
            main = method, xlab = paste(method, "1"), ylab = paste(method,
            "2"), xlim = c(min(coordinates$x) * 1.5, max(coordinates$y) *
            1.1))
        abline(h = 0)
        abline(v = 0)
        symbols(x = coordinates$x, coordinates$y, circles = rep(1,
            nrow(coordinates)), inches = 0.04, bg = cr[cluster.id[,
            1]], add = TRUE)
        if (any(cluster.id[, 1] != cluster.id[, 2])) {
            legend("topleft", legend = obj@cluster$names, fill = cr,
                cex = 0.7)
        }
    }
    obj
})

# Cell clustering
# ===========================================
if (!isGeneric("cellClustering")) {
    if (is.function("cellClustering"))
        fun <- cellClustering 
    else fun <- function(obj, ...) standardGeneric("cellClustering")
    setGeneric("cellClustering", fun)
}
#' Clustering
#'
#' Identifies the cell clusters, i.e. the cell subpopulations.
#' @name cellClustering
#' @aliases cellClustering,SCSRDataModel-method
#'
#' @param obj an object of class SCSRDataModel
#' @param n.cluster a number, an estimation of the ideal number of 
#' clusters is computed if equal to 0
#' @param n a number, the maximum to consider for an automatic 
#' determination of the ideal number of clusters
#' @param projection.method 'tsne'
#' @param method 'kmeans' or 'simlr'
#' @param markers a data frame of cell type signature genes
#' @param classification the method to use for annotation 
#' (by cluster or a library), 'none' to not do annotation
#' @param plot a logical
#' @param verbose a logical
#' @param most.variables a logical
#' @param pdf a logical
#' @param write a logical
#'
#' @details If the user knows the number of clusters present in her data set,
#' then `n.cluster` can be set and the estimation of the number of clusters is
#' skipped. `n` is the maximum number of clusters that the automatic estimation
#' of the number of clusters will consider. 
#' @details It is ignored if `n.cluster` is provided. `method` must be 
#' 'simlr' or 'kmeans' exclusively. If set to 'simlr', then the function uses
#' the **SIMLR()** function (**SIMLR** package) to perform clustering.
#' @details If set to 'kmeans' the function will perform a
#' dimensionality reduction by principal component analysis (PCA) followed by
#' K-means clustering and 2-dimensional projection by t-distributed stochastic
#' neighbor embedding (t-SNE). Regardless of the value of `method` ('simlr' or
#' 'kmeans'), in case `n.cluster` is not provided, then the function relies on
#' the **SIMLR_Estimate_Number_of_Clusters()** function to determine the number
#' of clusters, between 2 and `n`. 
#' @details If `plot` is TRUE, then the function displays
#' the t-SNE map with each cell colored according to the cluster it belongs to.
#' @details If `method` argument is 'simlr', then it further displays a heatmap
#' of thesimilarity matrix calculated by the **SIMLR()** function. 
#' @details If `pdf` is TRUE, then the function exports
#' the t-SNE plot in a pdf file in the *images* folder. The file is named
#' 't-SNE_map-X.pdf', where X is the `method` argument.
#' @details If `write` is TRUE, then the function writes two text files in the
#' *data* folder. The first one is called 'cluster-Y-X.txt', containing the
#' cluster vector assigning each cell of `data` to a cluster. The second one is
#' called 'tsne-Y-X.txt', containing the coordinates of each cell in the 2D
#' t-SNE projection. 'X' is the `method` argument anf 'Y' is the retained 
#' number of clusters. 
#' @details If `most.variables` is TRUE, then the function uses the most 
#' variable genes matrix counts if it exists in the object.
#' @details If the user want to do cluster annotation the classification method
#' and celltype marker table must be provided.
#' @details The ` markers` argument must be a table with cell type gene 
#' signatures, one cell type in each column. The column names are the names 
#' of the cell types.
#' @details The *markers.default* table provides an example of this format.
#' @details The`classification` argument must be one of 'by cluster', or 
#' 'library'. 'by cluster' will assign each cluster to a celltype by comparing 
#' the overall expression with markers given in ` markers`
#' 'library' will use the package *insert package name* to perform the 
#' labelling.
#' @details If `write` is TRUE and classification is not equal to 'none', 
#' then the function writes in a second text file a table containing 
#' probabilities of assignments of each cluster to a cell type for each cell
#' cluster. This cell type calling is performed as for the individual cells 
#' without thresholding but based on the cluster average transcriptome.
#'
#' @return A SCSRDataModel with cluster definition
#'
#' @export
#' 
#' @import SIMLR
#' @import Rtsne
#' @importFrom stats kmeans
#' @importFrom stats median
#' @importFrom stats model.matrix
#' @importFrom stats p.adjust
#' @importFrom stats phyper
#' @importFrom stats pnorm
#' @importFrom stats prcomp
#' @importFrom stats quantile
#' @importFrom stats sd
#' @importFrom stats var
#' @importFrom grDevices dev.off
#' @importFrom grDevices rainbow
#' @importFrom graphics abline
#' @importFrom graphics barplot
#' @importFrom graphics legend
#' @importFrom graphics symbols
#' @import pheatmap
#'
#' @examples
#' message('cellClustering')
#' message('--dataPrepare')
#' data <- matrix(runif(1000,0,1),nrow=50,ncol=20)
#' rownames(data) <- paste('gene',seq_len(50))
#' obj <- dataPrepare(data)
#' message('--cell Clustering')
#' obj <- cellClustering(obj)

setMethod("cellClustering", "SCSRDataModel", function(obj, n.cluster = 0,
    n = 10, most.variables = TRUE, projection.method = "tsne",
    method = c("simlr", "kmeans"), classification = c("none",
        "by cluster", "singler"), markers = markers_default,
    verbose = TRUE, plot = TRUE, pdf = FALSE, write = FALSE) {


    data <- obj@ncounts$matrix

    if (!is.null(obj@ncounts$matrix.mv) & most.variables) {
        data <- obj@ncounts$matrix.mv
        message("Matrix of most variable genes used. To use the whole matrix  
            set most.variables parameter to FALSE.")
    }

    if (!is.null(obj@cluster$id)) {
        message("Cluster definition will be overwritten.")
    }

    if (!dir.exists("images") & pdf) {
        dir.create("images")
    }
    if (!dir.exists("data") & write) {
        dir.create("data")
    }
    if (n.cluster == 0) {
        message("Estimating the number of clusters")
        c <- SIMLR_Estimate_Number_of_Clusters(data, NUMC = 2:n)
        a <- data.matrix(as.numeric((c$K1 + c$K2)))
        rownames(a) <- c(2:n)
        n.cluster <- as.numeric(rownames(subset(a, a == min(a))))
        message("Estimated number of clusters = ", n.cluster)
    }
    method <- match.arg(method)

    final <- list()
    if (method == "simlr") {
        s <- SIMLR(data, c = n.cluster, no.dim = 2)
        Y <- s$y
        cluster <- Y$cluster
        final[[1]] <- cluster
        final[[2]] <- s$F
        final[[3]] <- Y$size
        final[[4]] <- log(s$S * 10^6 + 1)
        names(final) <- c("cluster", "t-SNE", "numbers", "similarity")
    }
    if (method == "kmeans") {
        data <- data[rowSums(data) > 0, ]
        pca <- prcomp(t(data), center = TRUE, scale. = TRUE)
        v <- (pca$sdev^2)/(sum(pca$sdev^2))
        n <- min(which(diff(v) > -10^-4))
        if (n == 0 | identical(n, integer(0))) {
            n <- round(ncol(pca$x)/2)
        }
        pca <- pca$x[, seq_len(n)]
        tsne <- Rtsne(pca)
        km <- kmeans(pca, n.cluster)
        cluster <- km$cluster
        final[[1]] <- cluster
        final[[2]] <- tsne$Y
        final[[3]] <- km$size
        names(final) <- c("cluster", "t-SNE", "numbers")
    }

    classification <- match.arg(classification)
    cr <- rainbow(max(cluster))
    if (plot & classification == "none") {
        if (method == "simlr") {
            pheatmap(log(s$S * 10^6 + 1))
        }
        plot(x = final[[2]][, 1], y = final[[2]][, 2], type = "n",
            main = "t-SNE Map", xlab = "t-SNE1", ylab = "t-SNE2",
            xlim = c(min(final[[2]][, 1]) * 1.5, max(final[[2]][,
                1]) * 1.1))
        abline(h = 0)
        abline(v = 0)
        symbols(x = final[[2]][, 1], y = final[[2]][, 2], circles = rep(1,
            nrow(final[[2]])), inches = 0.04, bg = cr[cluster],
            add = TRUE)
        legend("topleft", legend = paste("cluster", seq_len(max(cluster))),
            fill = cr, cex = 0.7)
    }
    if (pdf) {
        pdf(paste("./images/t-SNE_map-", method, ".pdf", sep = ""))
        plot(x = final[[2]][, 1], y = final[[2]][, 2], type = "n",
            main = "t-SNE Map", xlab = "t-SNE1", ylab = "t-SNE2",
            xlim = c(min(final[[2]][, 1]) * 1.5, max(final[[2]][,
                1]) * 1.1))
        abline(h = 0)
        abline(v = 0)
        symbols(x = final[[2]][, 1], y = final[[2]][, 2], circles = rep(1,
            nrow(final[[2]])), inches = 0.04, bg = cr[cluster],
            add = TRUE)
        legend("topleft", legend = paste("cluster", seq_len(max(cluster))),
            fill = cr, cex = 0.7)
        dev.off()
    }
    if (write) {
        fwrite(data.frame(final[[1]]), paste("./data/cluster-",
            n.cluster, "-", method, ".txt", sep = ""), sep = "\t")
        fwrite(data.frame(final[[2]]), paste("./data/tsne-",
            n.cluster, "-", method, ".txt", sep = ""), sep = "\t")
    }
    message(n.cluster, " clusters detected")
    for (i in seq_len(n.cluster)) {
        message("cluster ", i, " -> ", final[[3]][i], " cells")
    }

    obj@cluster$method <- method
    obj@cell.representation$coordinates <- final[["t-SNE"]]
    obj@cell.representation$method <- projection.method
    obj@cluster$id <- final$cluster
    names <- sort(unique(final$cluster), decreasing = FALSE)
    names(names) <- names
    obj@cluster$names <- names

    if (classification != "none") {
        message("Annotation of calculated clusters.")
        obj <- .clusterClassifying(obj, markers, classification,
            most.variables, plot, write, verbose)
    }
    obj
})

#' Cluster classifier
#'
#' Classifies each cluster using cell type specific markers.
#'
#' @param obj an object of type SCSRDataModel
#' @param markers a data frame of cell type signature genes
#' @param method the method to use for annotation (by cluster or a library)
#' @param most.variables a logical
#' @param plot a logical (if TRUE, then plots the number of cells
#' attributed to one cell type, see below)
#' @param write a logical
#' @param verbose a logical
#' @keywords internal
#' 
#' @importFrom celldex HumanPrimaryCellAtlasData
#' @importFrom SingleR SingleR
#'
#' @return A SCSRDataModel with cluster definition

.clusterClassifying <- function(obj, markers, method, most.variables,
    plot, write, verbose) {

    coord <- obj@cell.representation$coordinates
    coord.method <- obj@cell.representation$method
    data <- obj@ncounts$matrix

    if (initialOrganism(obj) != "hsapiens")
        genes <- unlist(obj@ncounts$initial.orthologs) 
    else genes <- rownames(data)

    if (!is.null(obj@ncounts$matrix.mv) & most.variables) {
        data <- obj@ncounts$matrix.mv
        if (initialOrganism(obj) != "hsapiens")
            genes <- unlist(obj@ncounts$initial.orthologs.mv) 
        else genes <- rownames(data)
    }

    if (!dir.exists("cell-classification") & write) {
        dir.create("cell-classification")
    }

    if (is.null(obj@cluster$id)) {
        message("Please assign each cell to a cluster before doing celltype  
            annotation (for example by using cellClustering function).")
    } else {
        # Method to assign each cell to a celltype through
        # markers table, then conclude for each cluster by
        # assigning it to the celltype in majority
        if (method == "by cluster") {
            rownames(data) <- genes
            res <- .classifyByCluster(data, obj@cluster$id, obj@cluster$names,
                markers, write, verbose)
            c.names <- unlist(res[["class"]])
            names(c.names) <- sort(unique(obj@cluster$id), decreasing = FALSE)
            obj@cluster$names <- c.names
            obj@cluster$markers <- markers
        } else {
            ### Add option to do the celltype annotation
            ### with a specific library

            hpca.se <- celldex::HumanPrimaryCellAtlasData()
            pred.hesc <- SingleR(test = data, ref = hpca.se,
                assay.type.test = 1, labels = hpca.se$label.main,
                de.method = "wilcox", clusters = obj@cluster$id)
            names(pred.hesc$labels) <- rownames(pred.hesc)
            obj@cluster$names <- pred.hesc$labels
            obj@cluster$markers <- "HumanPrimaryCellAtlasData"

        }

        if (length(obj@cluster$names) != max(obj@cluster$id) |
            grepl("/", paste(obj@cluster$names, collapse = ""))) {
            stop("The length of obj@cluster$names must be equal to the number 
                of clusters and must contain no duplicates. The cluster names  
                must not include special characters")
        } else if (sum(duplicated(obj@cluster$names)) > 0) {
            message("Several clusters were assigned to the same celltype: 
                grouping the clusters as one.")
            rep <- data.frame(id = as.numeric(names(obj@cluster$names)),
                name = obj@cluster$names)
            uniq <- data.frame(new.id = seq(1, 
                length(unique(obj@cluster$names)), 1), 
                name = unique(obj@cluster$names))
            rep <- merge(rep, uniq, by.x = "name")

            names <- as.character(uniq$name)
            names(names) <- uniq$new.id
            obj@cluster$names <- names

            num <- data.frame(id = obj@cluster$id, index = seq(1,
                length(obj@cluster$id), 1))
            num <- merge(num, rep, by.x = "id")
            num <- num[order(num$index), ]
            obj@cluster$id <- num$new.id
        }


        obj@cluster$method <- method

        cr <- rainbow(max(obj@cluster$id))

        if (!is.null(coord)) {
            plot(x = coord[, 1], y = coord[, 2], type = "n",
                main = paste(coord.method, "by cluster"), 
                xlab = paste(coord.method, "1"), 
                ylab = paste(coord.method, "2"), 
                xlim = c(min(coord[, 1]) * 2, max(coord[, 1]) * 1.05))
            abline(h = 0)
            abline(v = 0)
            symbols(x = coord[, 1], y = coord[, 2], circles = rep(1,
                nrow(coord)), inches = 0.04, bg = cr[obj@cluster$id],
                add = TRUE)
            legend("topleft", legend = obj@cluster$names, fill = cr,
                cex = 0.75)
        }
    }
    obj
}


#' Cluster classification
#'
#' Analysis of the clusters composition by a marker based approach.
#'
#' @param data A counts matrix
#' @param cluster A vector of cluster ids
#' @param c.names A vector with unique cluster names
#' @param markers a table of cell type signature genes
#' @param write a logical
#' @param verbose a logical
#' @keywords internal

.classifyByCluster <- function(data, cluster, c.names, markers,
    write, verbose) {

    genes <- rownames(data)

    if (!dir.exists("cluster-analysis") & write) {
        dir.create("cluster-analysis")
    }
    if (is.null(c.names)) {
        c.names <- paste("cluster", seq_len(max(cluster)))
    }
    if (min(cluster) != 1) {
        cluster <- cluster + 1 - min(cluster)
    }

    n.cluster <- max(cluster)

    class <- list()
    length(class) <- n.cluster
    z <- c(seq_len(n.cluster))
    final.2 <- NULL
    types <- NULL

    final <- matrix(0, nrow = n.cluster, ncol = ncol(markers))
    colnames(final) <- colnames(markers)
    rownames(final) <- c.names
    for (i in z) {
        tmp <- NULL
        for (j in seq_len(ncol(markers))) {
            m.genes <- markers[, j][markers[, j] %in% genes]
            tmp <- c(tmp, sum(rowSums(data.frame(data[m.genes,
                cluster == i]))/sum(cluster == i))/length(m.genes))
            tmp[is.na(tmp)] <- 0
        }
        final[i, ] <- as.numeric(tmp)
    }
    final <- final/rowSums(final, na.rm = TRUE)
    final.2 <- final
    m <- apply(final.2, 1, function(x) (max(x) - min(x))/2)
    for (i in z) {
        final.2[i, ][final.2[i, ] < m[i]] <- 0
        c.type <- names(which(final.2[i, ] != 0))
        if (length(c.type) == 1) {
            if (verbose) {
                message(c.names[i], ": ", c.type)
            }
            class[[i]] <- c.type
        }
        if (length(c.type) == 2) {
            if (verbose) {
                message(c.names[i], ": ", c.type[1], " or ",
                    c.type[2])
            }
            class[[i]] <- paste(c.type[1], " or ", c.type[2],
                sep = "")
        }
        if (length(c.type) > 2) {
            if (verbose) {
                message(c.names[i], ": Unconclusive (more than 2 cell types
                    attributed to the cluster)")
            }
            class[[i]] <- "Undefined"
        }
    }
    names(class) <- c.names
    pheatmap(t(final.2), cluster_cols = FALSE)
    types <- do.call(rbind, class)
    if (write) {
        fwrite(data.frame(final), paste("./cluster-analysis/cluster_types.txt",
            sep = ""), sep = "\t")
    }
    res <- list(class, types, final.2)
    names(res) <- c("class", "types", "matrix")

    res
}

# Cell classification
# ===========================================
if (!isGeneric("cellClassifying")) {
    if (is.function("cellClassifying"))
        fun <- cellClassifying 
    else fun <- function(obj, ...) standardGeneric("cellClassifying")
    setGeneric("cellClassifying", fun)
}

#' Cell classifier
#'
#' Classifies each cell using cell type specific markers.
#' @name cellClassifying
#' @aliases cellClassifying,SCSRDataModel-method
#'
#' @param obj an object of type SCSRDataModel
#' @param markers a data frame of cell type signature genes
#' @param projection.method 'tsne'
#' @param plot a logical (if TRUE, then plots the number of cells
#' attributed to one cell type, see below)
#' @param most.variables a logical
#' @param write a logical
#' @param verbose a logical
#'
#' @details The ` markers` argument must be a table with cell type gene
#' signatures, one cell type in each column. The column names are the 
#' names of the cell types.
#' @details The *markers.default* table provides an example of this format.
#' @details If no projection is provided in the object, then the function 
#' will not display the cells. 
#' @details If ` plot` is TRUE, then the function plots the number of cells
#' attributed to a single cell type as a function of the threshold applied to
#' the normalized gene signature average.
#' @details If `most.variables` is TRUE, then the function uses the most
#' variable genes matrix counts if it exists in the object.
#' @details If ` write` is TRUE, then the function writes four different
#' text files. (1) The 'raw classification matrix' provides the normalized
#' average gene signature for each cell type in each individual cell, a number
#' between 0 and 1. This matrix has one row per cell type and one column per
#' cell, and the sum per column is 1. Row names are the cell type names (column
#' names of the markers table) and the column names are the individual cell
#' identifiers (column names of `data`). (2) The 'thresholded classification
#' matrix', which is obtained by eliminating all the values of the 'raw
#' classification matrix' that are below a threshold a\*. In practice, a\* is
#' automatically determined by the function to maximize the number of cells
#' that are assigned to a single cell type and all the cells (columns) assigned
#' to 0 or >1 cell types are discarded. The number of cells assigned to a
#' single type depending on a\* can be plotted by using the parameter
#' `plot=TRUE`. (3) A cluster vector assigning each cell to a cell type. Note
#' that a supplementary, virtual cluster is created to collect all the cells
#' assigned to 0 or >1 types. This virtual cluster is named 'undefined'.
#' (4) A table associating each cell type to a cluster number 
#' in the cluster vector.
#'
#' @return A SCSRDataModel with cluster definition
#'
#' @export
#'
#' @examples
#' message('cellClassifying')
#' message('--dataPrepare')
#' data <- matrix(runif(1000,0,1),nrow=50,ncol=20)
#' rownames(data) <- paste('gene',seq_len(50))
#' obj <- dataPrepare(data)
#' message('--cell Classifying')
#' markers <- matrix(paste('gene',seq_len(50)),ncol=5,nrow=10)
#' colnames(markers) <- paste('type',seq_len(5))
#' obj <- cellClassifying(obj, markers = markers)

setMethod("cellClassifying", "SCSRDataModel", function(obj, 
    markers = markers_default, projection.method = "tsne", 
    most.variables = TRUE, plot = TRUE, write = FALSE, verbose = TRUE) {

    coord <- obj@cell.representation$coordinates
    coord.method <- obj@cell.representation$method
    data <- obj@ncounts$matrix


    if (!is.null(obj@cluster$id)) {
        message("Cluster definition will be overwritten.")
    }
    obj@cluster$id <- NULL
    obj@cluster$names <- NULL

    if (initialOrganism(obj) != "hsapiens")
        genes <- unlist(obj@ncounts$initial.orthologs) 
    else genes <- rownames(data)

    if (!is.null(obj@ncounts$matrix.mv) & most.variables) {
        message("Matrix of most variable genes used. To use the whole matrix
            set most.variables parameter to FALSE.")
        data <- obj@ncounts$matrix.mv
        if (initialOrganism(obj) != "hsapiens")
            genes <- unlist(obj@ncounts$initial.orthologs.mv) 
        else genes <- rownames(data)
    }

    if (!dir.exists("cell-classification") & write) {
        dir.create("cell-classification")
    }

    if (is.null(coord) & projection.method == "tsne") {
        message("Calculation of 2D projection.")
        data2 <- data[rowSums(data) > 0, ]
        pca <- prcomp(t(data2), center = TRUE, scale. = TRUE)
        v <- (pca$sdev^2)/(sum(pca$sdev^2))
        n <- min(which(diff(v) > -10^-4))
        if (n == 0 | identical(n, integer(0))) {
            n <- round(ncol(pca$x)/2)
        }
        if (n != Inf){
            pca <- pca$x[, seq_len(n)]
            tsne <- Rtsne(pca)
            coord <- tsne$Y
            coord.method <- projection.method
        }
        obj@cell.representation$coordinates <- coord
        obj@cell.representation$method <- coord.method
    }

    # Method to assign each cell to a celltype through
    # markers table, then conclude for each cluster by
    # assigning it to the celltype in majority
    rownames(data) <- genes

    res <- .classifyByCell(data, coord, coord.method, markers,
        plot, write, verbose)

    obj@cluster$id <- res[["cluster"]]
    c.names <- res[["c.names"]]
    names(c.names) <- seq_len(max(obj@cluster$id))
    obj@cluster$names <- c.names


    obj@cluster$method <- "by cell"


    obj
})


#' Cell classifier
#'
#' Classifies each cellga using cell type specific markers.
#'
#' @param data    A counts matrix.
#' @param coordinates    A dataframe with projection coordinates.
#' @param projection    A character string with the projection method.
#' @param markers    A dataframe with the markers of celltypes.
#' @param plot    A logical.
#' @param write    A logical.
#' @param verbose    A logical.
#' @keywords internal
#' 
#' @return Return a list with each cell assignment cluster and
#' other informations.  

.classifyByCell <- function(data, coordinates, projection, markers,
    plot, write, verbose) {
    
    n.types <- ncol(markers)
    genes <- rownames(data)

    tmp <- data[as.character(unlist(markers))[as.character(unlist(markers)) 
        %in% genes], ]

    if (is.null(dim(tmp))) {
        stop("Not enough markers genes to pursue the cell classification")
    }
    final <- matrix(0, ncol = ncol(tmp), nrow = (ncol(markers)))

    for (i in seq_len(n.types)) {
        m.genes <- markers[, i][markers[, i] %in% genes]
        if (!is.null(dim(tmp[m.genes, ])))
            final[i, ] <- colSums(tmp[m.genes, ])/length(m.genes) 
        else final[i, ] <- NA
    }
    final[is.na(final)] <- 0
    final[, colSums(final) != 0] <- apply(final[, colSums(final) !=
        0], 2, function(x) x/sum(x))
    # final <- round(final*1000)/10
    rownames(final) <- colnames(markers)
    colnames(final) <- colnames(data)
    l <- matrix(0, 101, 2)
    q <- 0
    # Counts the number of cells for which only one
    # celltype is prevalent at a certain threshold
    for (n in seq(0.01, 1, 0.01)) {
        q <- q + 1
        f <- final
        f[f < n] <- 0
        l[q + 1, 1] <- n
        l[q + 1, 2] <- sum(apply(f, 2, function(x) sum(x == 0) ==
            (n.types - 1)))
    }
    # Select the max threshold for which the highest number
    # of cells is assigned to only one celltype
    seuil <- max(l[which(l[, 2] == max(l[, 2])), 1])
    m <- matrix(0, (n.types + 1), 2)
    for (n in 0:n.types) {
        f <- final
        f[f < seuil] <- 0
        f <- matrix(f[, apply(f, 2, function(x) sum(x == 0) ==
            n)], nrow = n.types)
        m[n + 1, 1] <- n
        m[n + 1, 2] <- sum(apply(f, 2, function(x) sum(x == 0) ==
            n))
    }
    m <- matrix(m[m[, 2] != 0, ], ncol = 2)

    res <- list()
    final.s <- final
    final.s[final.s < seuil] <- 0
    # Select columns (ie. cells) for which only one
    # celltype is assigned with the selected threshold
    res[[1]] <- final.s[, apply(final.s, 2, function(x) sum(x ==
        0) == (n.types - 1))]
    res[[3]] <- final
    final[!apply(final, 2, function(x) x == max(x))] <- 0
    res[[2]] <- final
    h <- pheatmap::pheatmap(res[[1]], cluster_rows = TRUE, cluster_cols = TRUE,
        show_colnames = FALSE)

    # Number of celltypes assigned directly to at least one
    # cell
    nn <- seq_len(length(rownames(res[[1]])[rowSums(res[[1]]) !=
        0]))
    names(nn) <- rownames(res[[1]])[rowSums(res[[1]]) != 0]
    # Assign a cluster to each cell for which only one
    # celltype could be assigned
    cluster <- unlist(lapply(apply(final.s, 2, function(x) names(x[x >
        0])), function(x) if (length(x) == 1) {
            return(nn[as.character(x)])
        } else {
            return(0)
        }))
    cluster[cluster == 0] <- max(cluster) + 1
    res[[4]] <- cluster

    if (sum(m[, 2] != 0) == 1) {
        n <- names(nn)
    } else {
        n <- c(names(nn), "Undefined cells")
    }
    d <- data.frame(Cell_type = c(rownames(res[[1]]), "Undefined"),
        Number_of_cells = c(apply(res[[1]], 1, function(x) sum(x !=
            0)), ncol(res[[3]]) - ncol(res[[1]])))
    d <- d[d$Number_of_cells != 0, ]
    rownames(d) <- paste("cluster", seq_len(nrow(d)))

    res[[5]] <- d$Cell_type

    if (verbose) {
        message("A threshold of ", round(seuil * 1000)/10, "% maximizes the 
            number of cells assigned to one cell type")
        for (i in (nrow(m)):1) {
            message(m[i, 2], " cell(s) or ", round(m[i, 2] *
                1000/ncol(data))/10, "% identified to ", (n.types) -
                m[i, 1], " cell type(s)")
        }
        message(" ")
    }
    if (sum(d$Number_of_cells == 1) > 0) {
        d$NewCell_type <- d$Cell_type
        d$NewCell_type[d$Number_of_cells == 1] <- "Undefined"
        d$Id_ini <- seq(1, nrow(d), 1)

        if ("Undefined" %in% d$NewCell_type) {
            name_new <- c(unique(d$NewCell_type[d$NewCell_type !=
                "Undefined"]), "Undefined")
        } else {
            name_new <- unique(d$NewCell_type)
        }
        new <- data.frame(NewCell_type = name_new, NewId = seq(1,
            length(unique(d$NewCell_type)), 1))

        d <- merge(d, new, by.x = "NewCell_type", sort = FALSE)
        d <- d[order(d$Id_ini), ]

        cluster_df <- data.frame(Id_ini = cluster)
        cluster_df$id <- seq(1, nrow(cluster_df), 1)
        cluster_df <- merge(cluster_df, d, by.x = "Id_ini", sort = FALSE)
        cluster_df <- cluster_df[order(cluster_df$id), ]
        cluster <- cluster_df$NewId
        c.names <- new$NewCell_type
        names(c.names) <- new$NewId

        message(max(d$Id_ini) - max(d$NewId), " clusters contained only 1 cell
            and are reassigned to 'Undefined'")

        d[d$Cell_type == "Undefined", "Number_of_cells"] <- 
            sum(d$Number_of_cells[d$NewCell_type ==  "Undefined"])

        d <- d[d$Cell_type %in% c.names, c("Cell_type", "Number_of_cells")]

        rownames(d) <- paste("cluster", seq_len(nrow(d)))
        res[[4]] <- cluster
        res[[5]] <- c.names

        n <- c.names
    }


    cr <- c(rainbow(max(cluster) - 1), "gray")


    names(res) <- c("tresh_mat", "max_mat", "raw_mat", "cluster",
        "c.names")

    if (verbose)
        message(paste0(capture.output(d), collapse = "\n"))

    if (!is.null(coordinates)) {
        plot(x = coordinates[, 1], y = coordinates[, 2], type = "n",
            main = "t-SNE Map Cell by cell", xlab = "t-SNE1",
            ylab = "t-SNE2", xlim = c(min(coordinates[, 1]) *
                2, max(coordinates[, 1]) * 1.05))
        abline(h = 0)
        abline(v = 0)
        symbols(x = coordinates[, 1], y = coordinates[, 2], circles = rep(1,
            nrow(coordinates)), inches = 0.04, bg = cr[cluster],
            add = TRUE)
        legend("topleft", legend = n, fill = cr, cex = 0.75)
    }
    if (plot) {
        plot(l, xlab = "Threshold (%)", ylab = "Number of cells",
            main = "Attribution of one cell to one cell type",
            type = "l")
        abline(v = l[l[, 1] == seuil, 1], col = "red", lty = 2)
        symbols(l[l[, 1] == seuil, 1], l[l[, 1] == seuil, 2],
            circles = 1, inches = 0.05, bg = "red", fg = "red",
            add = TRUE)
        text(x = l[l[, 1] == seuil, 1], y = l[l[, 1] == seuil,
            2] - 0.07 * (max(l)), labels = paste(seuil, "%"),
            offset = 20)
    }
    if (write) {
        fwrite(data.frame(cbind(rownames(res[[1]]), res[[1]])),
            "cell-classification/threshold_class_matrix.txt",
            sep = "\t")
        fwrite(data.frame(cbind(rownames(res[[3]]), res[[3]])),
            "cell-classification/raw_class_matrix.txt", sep = "\t")
        nn <- cbind(names(nn), nn)
        colnames(nn) <- c("cluster.name", "cluster.number")
        fwrite(data.frame(nn), "cell-classification/cluster2name.txt",
            sep = "\t")
        fwrite(data.frame(cluster), "cell-classification/cluster.txt",
            sep = "\t")
        c.names <- res[[5]]
    }
    res
}

if (!isGeneric("dgeCluster")) {
    if (is.function("dgeCluster"))
        fun <- dgeCluster 
    else fun <- function(obj, ...) standardGeneric("dgeCluster")
    setGeneric("dgeCluster", fun)
}

#' Differential Gene Expression
#'
#' Analysis of the differentially expressed genes in the clusters
#' @name dgeCluster
#' @aliases dgeCluster,SCSRDataModel-method
#'
#' @param obj A SCSRDataModel object
#' @param s.pval A pvalue threshold
#' @param most.variables a logical
#' @param write A logical
#' @param verbose A logical
#' @param plot A logical
#' 
#' @details Computes the diferential gene expression
#' between the clusters using **edgeR**)
#' @details If `most.variables` is TRUE, then the function uses the 
#' most variable genes matrix counts if it exists in the object.
#' @details If `write` is TRUE, then the
#' function writes a text file named 'table_dge_X.txt', where X is the 
#' cluster name, that contains the list of differentially expressed genes. 
#'
#' @return A SCSRDataModel with cluster definition
#'
#' @export
#' 
#' @import edgeR
#' @import limma
#' @import RColorBrewer
#'
#' @examples
#' message('dgeCluster')
#' message('--dataPrepare')
#' data <- matrix(runif(1000,0,1),nrow=50,ncol=20)
#' rownames(data) <- paste('gene',seq_len(50))
#' obj <- dataPrepare(data)
#' message('--add Clustering')
#' obj <- addClustering(obj,cluster.id = sample(1:5,size = 20,replace = TRUE))
#' if (FALSE) obj <- dgeCluster(obj)

setMethod("dgeCluster", "SCSRDataModel", function(obj, s.pval = 10^-2,
    most.variables = TRUE, write = FALSE, verbose = TRUE, plot = TRUE) {


    if (!dir.exists("cluster-analysis") & write) {
        dir.create("cluster-analysis")
    }

    data <- obj@ncounts$matrix

    if (initialOrganism(obj) != "hsapiens")
        genes <- unlist(obj@ncounts$initial.orthologs) 
    else genes <- rownames(data)

    if (!is.null(obj@ncounts$matrix.mv) & most.variables) {
        message("Matrix of most variable genes used. To use the whole matrix 
        set most.variables parameter to FALSE.")
        data <- obj@ncounts$matrix.mv
        if (initialOrganism(obj) != "hsapiens")
            genes <- unlist(obj@ncounts$initial.orthologs.mv) 
        else genes <- rownames(data)
    }
    cluster <- obj@cluster$id
    c.names <- obj@cluster$names

    if (!is.null(obj@ncounts$matrix.mv) & most.variables) {
        message("Matrix of most variable genes used. To use the whole matrix   
            set most.variables parameter to FALSE.")
        data <- obj@ncounts$matrix.mv
    }

    genes <- rownames(data)
    diff.genes <- list()
    z <- seq(length(c.names))
    k <- 0
    v <- NULL

    dge <- DGEList(data, genes = genes)
    dge <- calcNormFactors(dge)

    if (verbose) {
        message("edgeR differential gene expression (dge) processing:")
    }

    for (i in z) {
        if (verbose) {
            message("Looking for differentially expressed genes in ",
                c.names[i])
        }
        cl <- as.matrix(cluster)
        cl[cluster == i, ] <- 1
        cl[cluster != i, ] <- 2
        c <- factor(cl)
        design <- model.matrix(~0 + c)
        cm <- makeContrasts(c1 - c2, levels = design)
        y <- estimateDisp(dge, design, robust = TRUE)
        fit.y <- glmFit(y, design)
        lrt <- glmLRT(fit.y, contrast = cm)
        sel.r <- topTags(lrt, p.value = s.pval, adjust.method = "BH",
            n = nrow(y$counts))
        if (length(sel.r) == 0) {
            if (verbose == TRUE) {
                message("No differentially expressed genes in ",
                    c.names[i])
            }
            v <- c(v, i)
        } else {
            k <- k + 1
            resu <- sel.r$table
            resu <- resu[order(resu$logFC, decreasing = TRUE), ]
            diff.genes[[k]] <- resu
            if (write == TRUE) {
                fwrite(data.frame(resu), paste("./cluster-analysis/table_dge_",
                    c.names[i], ".txt", sep = ""), sep = "\t")
            }
        }
    }
    if (is.null(v)) {
        names(diff.genes) <- c.names
    } else {
        names(diff.genes) <- c.names[-v]
    }

    obj@dge.cluster$genes <- diff.genes
    obj@dge.cluster$param <- list(s.pval = s.pval)

    if (plot) {
        message("Plotting the 50 genes with highest log 
            fold change for each cluster.")
        clustersAnnot <- data.frame(id = cluster, index = seq_len(
            length(cluster)))
        annot <- data.frame(id = names(c.names), name = c.names)
        clustersAnnot <- merge(clustersAnnot, annot, by.x = "id",
            order = FALSE)
        clustersAnnot <- clustersAnnot[order(clustersAnnot$id), ]
        clustersAnnot <- data.frame(cluster = clustersAnnot$name)

        if (length(diff.genes)>0){
            pdata <- data[unlist(lapply(diff.genes, function(x) 
                x[order(x$logFC, decreasing = TRUE), ]
                [seq_len(min(nrow(x), 50)), "genes"])), 
                order(cluster)]
            rownames(clustersAnnot) <- colnames(pdata)
            if (length(c.names) <= 2) {
                col <- brewer.pal(3, "Paired")
                col <- col[seq_len(length(c.names))]
            } else col <- brewer.pal(length(c.names), "Paired")
                names(col) <- c.names
                color <- list(cluster = col)
                pheatmap(pdata, cluster_rows = FALSE, cluster_cols = FALSE,
                    show_rownames = FALSE, show_colnames = FALSE, 
                    main = "Main DE genes per cluster",
                    annotation_col = clustersAnnot, 
                    annotation_colors = color)
        }else{
            message("No differentially expressed genes in any cluster.")
        }
    }
    obj
})



if (!isGeneric("cellSignaling")) {
    if (is.function("cellSignaling"))
        fun <- cellSignaling 
    else fun <- function(obj, ...) standardGeneric("cellSignaling")
    setGeneric("cellSignaling", fun)
}

#' Cell Signaling
#'
#' Computes 'autocrine' or 'paracrine' interactions between cell
#' clusters.
#' @name cellSignaling
#' @aliases cellSignaling,SCSRDataModel-method
#'
#' @param obj an object of type SCSRDataModel
#' @param int.type 'autocrine' or 'paracrine'
#' @param s.score LRscore threshold
#' @param logFC a number, the log fold-change threshold for differentially
#' expressed genes
#' @param tol a tolerance parameter for balancing 'autocrine|paracrine'
#' interactions to the 'autocrine' or 'paracrine' group
#' @param most.variables a logical
#' @param addLR NULL or a dataframe with 2 columns
#' @param switchDB a dataframe with 2 columns
#' @param write a logical
#' @param verbose a logical
#'
#' @details `int.type` must be equal to 'paracrine' or 'autocrine' exclusively.
#' The 'paracrine' option looks for ligands expressed in cluster A and their
#' associated receptors according to LR*db* that are expressed in any other
#' cluster but A. These interactions are labelled 'paracrine'. The interactions
#' that involve a ligand and a receptor, both differentially expressed in their
#' respective cell clusters according to the **edgeR** analysis performed by 
#' the **cluster_analysis()** function, are labelled 'specific'. The 
#' 'autocrine' option searches for ligands expressed in cell cluster A and
#' their associated receptors also expressed in A. These interactions are
#' labelled 'autocrine'. Additionally, it searches for those associated
#' receptors in the other cell clusters (not A) to cover the part of the
#' signaling that is 'autocrine' and 'paracrine' simultaneously. These
#' interactions are labelled 'autocrine/paracrine'.
#' @details
#' The `tol` argument allows the user to tolerate a fraction of the cells in
#' cluster A to express the receptors in case `int.type='paracrine'`, that is
#' to call interactions that are dominantly paracrine though not exclusively.
#' Conversely, it allows the user to reject interactions involving receptors
#' that would be expressed by a small fraction of cluster A cells in case
#' `int.type='autocrine'`. By construction theassociation of these two 
#' options covers all the possible interactions and increasing the `tol`
#' argument allows the user to move interactions from 'autocrine' to 
#' 'paracrine'.
#' @details
#' `s.score` is the threshold on the LRscore. The value must lie in the [0;1]
#' interval, default is 0.5 to ensure confident ligand-receptor pair
#' identifications (see our publication). Lower values increase the number of
#' putative interactions while increasing the false positives. Higher values
#' do the opposite.
#' @details
#' `logFC` is a threshold applied to the log fold-change (logFC) computed for
#' each gene during the differential gene expression analysis. Its default 
#' value is log~2~(1.5) It further selects the differentially expressed genes
#' (>logFC) after the p-value threshold imposed in the function 
#' **cluster_analysis()** below.
#' @details 
#' If `most.variables` is TRUE, then the function uses the most variable genes
#' matrix counts if it exists in the object.
#' @details 
#' The `addLR` allows the user to add LR interaction that they want to study.
#' It must be a dataframe with one column 'ligand' and one column 'receptor'.
#' @details 
#' The `switchDB` allows the user to switch the LR database used to study.
#' interactions.
#' @details
#' If `write` is TRUE, then the function writes a text file that reports the
#' interactions in the *cell-signaling* folder. This file is a 4-column table:
#' ligands, receptors, interaction types ('paracrine', 'autocrine',
#' 'autocrine/paracrine' and 'specific'), and the associated LRscore.
#' @details
#' Remarks:
#' @details- This function can be used with any `data` table associated with
#' corresponding `genes` and `cluster` vectors, meaning that advanced users 
#' can perform their own data normalization and cell clustering upfront.
#' @details- In case the function **cluster_analysis()** was not executed,
#' this function would work but 'specific' interactions would not be 
#' annotated as such.
#'
#' @return A SCSRInference with LR interaction
#'
#' @export
#'
#' @examples
#' message('cellSignaling')
#' message('--dataPrepare')
#' data <- matrix(runif(1000,0,1),nrow=50,ncol=20)
#' rownames(data) <- paste('gene',seq_len(50))
#' obj <- dataPrepare(data)
#' message('--add Clustering')
#' obj <- addClustering(obj,cluster.id = sample(1:5,size = 20,replace = TRUE))
#' message('--cell Signaling')
#' obj.int <- cellSignaling(obj,int.type = 'paracrine')


setMethod("cellSignaling", "SCSRDataModel", function(obj, 
    int.type = c("paracrine", "autocrine"), s.score = 0.5, logFC = log2(1.5), 
    tol = 0, most.variables = TRUE, addLR = NULL, switchDB = NULL, 
    write = FALSE, verbose = TRUE) {

    data(LRdb)

    if (!is.null(switchDB)) {
        if (!is.data.frame(switchDB)) {
            message("The database you input must be a dataframe. Using LRdb..")
        } else {
            LRdb <- switchDB[, seq_len(2)]
        }
    }
    cluster <- obj@cluster$id
    c.names <- obj@cluster$names
    if (is.null(cluster)) {
        cluster <- obj@cell.classification$identical
        c.names <- obj@cell.classification$names
    }
    species <- obj@initial.organism
    data <- obj@ncounts$matrix
    param <- list(s.score = s.score, logFC = logFC, tol = tol)

    if (!is.null(obj@ncounts$matrix.mv) & most.variables) {
        message("Matrix of most variable genes used. To use the whole matrix   
            set most.variables parameter to FALSE.")
        data <- obj@ncounts$matrix.mv
    }

    if (!dir.exists("cell-signaling") & write) {
        dir.create("cell-signaling")
    }
    if (is.null(c.names)) {
        c.names <- paste("cluster", seq_len(max(cluster)))
    }
    if (min(cluster) != 1) {
        cluster <- cluster + 1 - min(cluster)
    }
    if (length(c.names) != max(cluster) | sum(duplicated(c.names)) >
        0 | grepl("/", paste(c.names, collapse = ""))) {
        stop("The length of c.names must be equal to the number of obj@cluster
            and must contain no duplicates. The obj@cluster names must not 
            include special characters")
    }
    int.type <- match.arg(int.type)

    genes <- rownames(data)

    z <- seq_len(max(cluster))
    if (species != "hsapiens") {
        if (!any(LRdb[, 1] %in% genes)) {
            ortho <- data.frame(HsapiensLig = rownames(obj@ncounts$matrix),
                speciesLig = obj@ncounts$initial.orthologs)
            lig <- data.frame(HsapiensLig = LRdb$ligand, 
                HsapiensRec = LRdb$receptor, id = seq(1, nrow(LRdb), 1))

            lig <- merge(lig, ortho, by.x = "HsapiensLig", order = FALSE)
            lig <- lig[order(lig$id), ]

            names(ortho) <- c("HsapiensRec", "speciesRec")
            lig <- merge(lig, ortho, by.x = "HsapiensRec", order = FALSE)

            lig <- lig[order(lig$id), ]
            LRdb <- lig[, c("speciesLig", "speciesRec")]
            names(LRdb) <- c("ligand", "receptor")
        }
    }


    if (!is.null(addLR)) {
        if (!is.data.frame(addLR) | ncol(addLR) < 2) {
            message("Please input the ligand-receptors interactions you want  
                to add as a two-column dataframe.")
        } else {
            if (!any(addLR[, 1] %in% genes)) {
                if (species != "hsapiens") {
                    ortho <- data.frame(Hsapiens = 
                        rownames(obj@ncounts$matrix),
                        species = obj@ncounts$initial.orthologs)
                    lig <- data.frame(HsapiensLig = addLR[, 1],
                        HsapiensRec = addLR[, 2], id = seq(1, nrow(addLR), 1))

                    lig <- merge(lig, ortho, by.x = "HsapiensLig",
                        order = FALSE)
                    lig <- lig[order(lig$id), ]

                    names(ortho) <- c("HsapiensRec", "speciesRec")
                    lig <- merge(lig, ortho, by.x = "HsapiensRec",
                        order = FALSE)

                    lig <- lig[order(lig$id), ]
                    addLR <- lig[, c("speciesLig", "speciesRec")]
                    names(addLR) <- c("ligand", "receptor")
                }
            }
            add <- cbind(addLR[, seq_len(2)], data.frame(matrix(NA,
                ncol = ncol(LRdb) - 2), nrow = nrow(addLR)))
            names(add) <- names(LRdb)
            LRdb <- rbind(LRdb, add)
            LRdb <- unique(LRdb)
        }
    }
    lig <- unique(LRdb$ligand)
    rec <- unique(LRdb$receptor)

    data <- data.frame(data)
    data <- data[rowSums(data) > 0, ]
    med <- sum(data)/(nrow(data) * ncol(data))

    out <- list()
    ## Autocrine -------------------
    if (int.type == "autocrine") {
        if (verbose) {
            message("Autocrine signaling: ")
        }
        k <- 0
        int <- NULL
        n.int <- NULL
        if (verbose) {
            message("Checking for cell/cell signaling:")
        }
        for (i in z) {
            if (sum(cluster == i) > 1) {
                tmp <- data[, cluster == i]
                tmp <- tmp[rowSums(tmp) > 0, ]
                if (sum(is.element(lig, rownames(tmp))) > 0) {
                    lig.tmp <- rownames(tmp)[is.element(rownames(tmp),
                        lig)]
                } else {
                    lig.tmp <- NULL
                }

                final.tmp <- LRdb[is.element(LRdb$ligand, lig.tmp),
                    seq_len(2)]
                final.tmp <- data.frame(final.tmp, 
                    as.character(rep("autocrine|paracrine",
                    sum(is.element(LRdb$ligand, lig.tmp)))))
                m.lig <- rowSums(tmp[unique(final.tmp[, 1]),
                    ])/sum(cluster == i)
                names(m.lig) <- unique(final.tmp[, 1])

                if (sum(is.element(rec, rownames(tmp))) > 0) {
                    rec.tmp <- rownames(tmp)[is.element(rownames(tmp[apply(tmp,
                        1, function(x) sum(x > 0)) > tol * ncol(tmp),
                        ]), rec)]
                } else {
                    rec.tmp <- NULL
                }

                for (j in z) {
                    if (sum(cluster == i) > 1) {
                        temp <- data[, cluster == j]
                        temp <- temp[rowSums(temp) > 0, ]
                        if (sum(is.element(rec, rownames(temp))) > 0) {
                            rec.temp <- rownames(temp)[is.element(
                                rownames(temp), rec)]
                        } else {
                            rec.temp <- NULL
                        }
                        rec.temp <- rec.temp[is.element(rec.temp,
                            rec.tmp)]
                        m.rec <- rowSums(data.frame(temp[rec.temp,
                            ]))/sum(cluster == j)
                        names(m.rec) <- rec.temp

                        final <- final.tmp[is.element(final.tmp$receptor,
                            rec.temp), ]
                        final <- cbind(final, .LRscore(m.lig[final$ligand],
                            m.rec[final$receptor], med))

                        colnames(final) <- c(c.names[i], c.names[j],
                            "interaction type", "LRscore")

                        if (i == j) {
                            final$`interaction type` <- "autocrine"
                        }

                        final <- final[final[, 4] > s.score, ]
                        final <- final[order(final[, 4], decreasing = TRUE), ]


                        if (nrow(final) > 0) {
                            k <- k + 1
                            out[[k]] <- final
                            if (verbose) {
                                message(nrow(final), " interactions from ",
                                c.names[i], " to ", c.names[j])
                            }
                            int <- c(int, paste(i, "-", j, sep = ""))
                            n.int <- c(n.int, paste(c.names[i], "-",
                                c.names[j], sep = ""))
                            gr <- graph_from_data_frame(final, 
                                directed = FALSE)
                            if (write) {
                                fwrite(data.frame(final), 
                                    paste("./cell-signaling/LR_interactions_",
                                    c.names[i], "-", c.names[j], "-", int.type,
                                    ".txt", sep = ""), sep = "\t")
                            }
                        }
                    }
                }
            }
        }
        if (k != 0) {
            names(out) <- n.int
        }
    }


    ## Paracrine -------------------
    if (int.type == "paracrine") {
        gene.list <- vector("list", max(cluster))
        for (i in z) {
            if (is.null(obj@dge.cluster$genes[[c.names[i]]]) &
                !file.exists(paste("./cluster-analysis/table_dge_",
                    c.names[i], ".txt", sep = ""))) {
                gene.list[[i]] <- "none"
                message("No such file as table_dge_", c.names[i],
                    ".txt in the cluster-analysis folder")
            } else {
                if (file.exists(paste("./cluster-analysis/table_dge_",
                    c.names[i], ".txt", sep = ""))) {
                    resu <- fread(paste("./cluster-analysis/table_dge_",
                        c.names[i], ".txt", sep = ""), data.table = FALSE)
                } else {
                    resu <- obj@dge.cluster$genes[[c.names[i]]]
                }
                gene.list[[i]] <- resu$genes[resu$logFC > logFC]

                gene.list[[i]] <- gene.list[[i]][!is.na(gene.list[[i]])]
            }
        }

        if (verbose) {
            message("Paracrine signaling: ")
        }
        k <- 0
        int <- NULL
        n.int <- NULL
        if (verbose) {
            message("Checking for signaling between cell types.")
        }
        for (i in z) {
            if (sum(cluster == i) > 1) {
                tmp <- data[, cluster == i]
                tmp <- tmp[rowSums(tmp) > 0, ]
                if (sum(is.element(lig, rownames(tmp))) > 0) {
                    lig.tmp <- rownames(tmp)[is.element(rownames(tmp),
                        lig)]
                } else {
                    lig.tmp <- NULL
                }

                final.tmp <- LRdb[is.element(LRdb$ligand, lig.tmp),
                    seq_len(2)]
                final.tmp <- data.frame(final.tmp, 
                    as.character(rep("paracrine",
                    sum(is.element(LRdb$ligand, lig.tmp)))), 
                    stringsAsFactors = FALSE)
                m.lig <- rowSums(tmp[unique(final.tmp[, 1]),
                    ])/sum(cluster == i)
                names(m.lig) <- unique(final.tmp[, 1])

                if (sum(is.element(rec, rownames(tmp))) > 0) {
                    rec.tmp <- rownames(tmp)[is.element(rownames(tmp[apply(tmp,
                        1, function(x) sum(x > 0)) > tol * ncol(tmp), ]), rec)]
                } else {
                    rec.tmp <- NULL
                }

                for (j in z[-i]) {
                    if (sum(cluster == j) > 1) {
                        temp <- data[, cluster == j]
                        temp <- temp[rowSums(temp) > 0, ]
                        if (sum(is.element(rec, rownames(temp))) > 0) {
                            rec.temp <- rownames(temp)[is.element(
                                rownames(temp), rec)]
                        } else {
                            rec.temp <- NULL
                        }
                        rec.temp <- rec.temp[!is.element(rec.temp,
                            rec.tmp)]
                        m.rec <- rowSums(data.frame(temp[rec.temp,
                            ]))/sum(cluster == j)
                        names(m.rec) <- rec.temp

                        final <- final.tmp[is.element(final.tmp$receptor,
                            rec.temp), ]
                        final <- cbind(final, .LRscore(m.lig[final$ligand],
                            m.rec[final$receptor], med))
                        exclus <- final$ligand %in% gene.list[[i]] &
                            final$receptor %in% gene.list[[j]]
                        if (sum(exclus) != 0) {
                            f.exclu <- final[exclus, ]
                            final <- final[!(final$ligand %in% 
                                gene.list[[i]] &
                                final$receptor %in% gene.list[[j]]), ]
                            f.exclu[, 3] <- "specific"
                            final <- rbind(f.exclu, final)
                        }

                        colnames(final) <- c(c.names[i], c.names[j],
                            "interaction type", "LRscore")
                        final <- final[final[, 4] > s.score, ]
                        final <- final[order(final[, 4], decreasing = TRUE), ]

                        if (nrow(final) > 0) {
                            k <- k + 1
                            out[[k]] <- final
                            if (verbose) {
                                message(nrow(final), " interactions from ",
                                c.names[i], " to ", c.names[j])
                            }
                            int <- c(int, paste(i, "-", j, sep = ""))
                            n.int <- c(n.int, paste(c.names[i], "-",
                                c.names[j], sep = ""))
                            gr <- graph_from_data_frame(final, 
                                directed = FALSE)
                            if (write) {
                                fwrite(data.frame(final), 
                                    paste("./cell-signaling/LR_interactions_",
                                    c.names[i], "-", c.names[j], "-", int.type,
                                    ".txt", sep = ""), sep = "\t")
                            }
                        } else {
                            if (verbose) {
                                message(nrow(final), 
                                    " No significant interaction found from ",
                                    c.names[i], " to ", c.names[j])
                            }
                        }
                    }
                }
            }
        }
        if (k != 0) {
            names(out) <- n.int
        }
    }
    if (length(out) > 0) {
        out_subpop <- data.frame(Subpop = names(out))
        out_subpop$Number <- unlist(lapply(out, nrow))
        LR <- foreach(int = seq_len(length(out)), .combine = "rbind") %do%
            {
                data.frame(Ligand = out[[int]][, 1], 
                    Receptor = out[[int]][, 2])
            }
        new("SCSRInference", LRinter = out, LRsubpop = out_subpop,
            ligands = LR$Ligand, receptors = LR$Receptor, param = param)
    } else {
        message("No interations found between the clusters. 
            Not returning an object.")
    }
})

#' Calculation of the LRscore
#'
#' @param l a value (or a vector) of the mean ligand expression in the
#' secreting cluster
#' @param r a value (or a vector) of the mean receptor expression in the
#' receiving cluster
#' @param s a value for scaling the score (usually the mean of the whole
#' read count table, the median or another similar value is possible),
#' must be over 
#' @keywords internal
#'
#' @return a value or a vector

.LRscore <- function(l, r, s) {
    L <- l^(1/2)
    R <- r^(1/2)
    S <- s
    sc <- L * R/(S + L * R)
    return(sc)
}
