#' @title inter network
#' @description Computes intercellular gene networks.
#'
#' @details
#' If `write` is TRUE, then the function writes four different files. A graphML
#' file in the *cell-signaling* folder for intercellular interactions between
#' each pair of clusters named "intercell_network_Z~1~-Z~2~.graphml", where
#' Z~1~ and Z~2~ are the *c.names* of the clusters. A graphML file in the
#' *cell-signaling* folder that contains a compilation of all the intercellular,
#' ligand-receptor interactions named "full-intercellular-network.graphml".
#' A text and a graphML file in the *networks* folder containing the intracellular
#' network for each cell cluster named "intracell_network_Z.txt" and
#' "intracell_network_Z.graphml", where Z is the *c.names* of the cluster.
#'
#' @param obj an object of type SCSRInference
#' @param dm an object of type SCSRDataModel
#' @param write a logical (if TRUE writes graphML and text files for the
#' interface and internal networks)
#' @param most.variables a logical
#' @param plot a logical
#' @param verbose a logical
#'
#' @return The function returns a list containing the tables of interaction
#' between two cell types and the table for the full network of all the cell
#' types.
#'
#' @export
#'
#' @importFrom igraph graph_from_data_frame set_vertex_attr write.graph V %>%
#' @import data.table
#'
#' @examples
#' print("dataPrepare")
#' data <- matrix(runif(1000,0,1),nrow=50,ncol=20)
#' rownames(data) <- paste("gene",seq_len(50))
#' obj <- dataPrepare(data)
#' print("cell Clustering")
#' obj <- cellClustering(obj)
#' print("cell Signaling")
#' obj.int <- cellSignaling(data,int.type = "paracrine")
#' net <- inter_network(obj.int, obj)

inter_network <- function(obj, dm, clusters = NULL, plot = TRUE,
          most.variables = TRUE, write = FALSE, verbose = TRUE){
  
  if (!is(dm, "SCSRDataModel")){
        stop("dm must be a SCSRDataModel object")
    }
    if (!is(obj, "SCSRInference")){
        stop("obj must be a SCSRInference object")
    }

  if (dir.exists("networks")==FALSE & write==TRUE){
    dir.create("networks")
  }

  c.names <- dm@cluster$names
  cluster <- dm@cluster$id
  data <- dm@ncounts$matrix
  genes <- rownames(dm@ncounts$matrix)
  signal <- obj@LRinter
  species <- dm@initial.organism
  if (species!='hsapiens'){
    rownames(data) <- dm@ncounts$initial.orthologs
  }

  if (!is.null(dm@ncounts$matrix.mv)&most.variables){
        cat("Matrix of most variable genes used. To use the whole matrix set most.variables 
            parameter to FALSE.\n")
        data <- dm@ncounts$matrix.mv
        genes <- rownames(dm@ncounts$matrix.mv)
        if (species!='hsapiens'){
          rownames(data) <- dm@ncounts$initial.orthologs.mv
      }
    }

  if (is.null(c.names)==TRUE){
    c.names <- paste("cluster",seq_len(max(cluster)))
  }
  if (min(cluster)!=1){
    cluster <- cluster + 1 - min(cluster)
  }
  if (length(c.names)!=max(cluster) | sum(duplicated(c.names))>0 |
      grepl("/",paste(c.names,collapse =""))){
    stop("The length of c.names must be equal to the number of clusters and
        must contain no duplicates. The cluster names must not include special
        characters")
  }
  cellint <- NULL
  rownames(data) <- genes
  interface <- list()

  clusterToStudy <- c.names
  if (!is.null(clusters)){
    if(all(clusters%in%c.names)){
      clusterToStudy <- clusters
      }else{
        cat("Please input clusters names that are described in the object. 
          Analyzing all clusters\n")
      }
  }
  

  # Interface networks ========================================================
  tmp <-  vector("list",length=length(signal))
  if (!is.null(signal)){
    for (i in seq_len(length(signal))){
      ci <-  signal[[i]]
      if(names(ci)[1]%in%clusterToStudy&names(ci)[2]%in%clusterToStudy){
        from <- names(ci)[1]
        to <- names(ci)[2]
        tmp[[i]] <-  data.frame(ligand=paste0(from,".",ci[[1]]),receptor=
                                paste0(to,".",ci[[2]]),ligand.name=ci[[1]],
                 receptor.name=ci[[2]],origin=from,destination=to,ci[,3:4],
                 stringsAsFactors=FALSE)
      }
    }
    cellint <- do.call("rbind",tmp)

    # Pair networks ----------------
    m <- 0
    n.int <- NULL
    for (i in seq_len(length(signal))){
      if(colnames(signal[[i]])[1]%in%clusterToStudy&
        colnames(signal[[i]])[2]%in%clusterToStudy){
        m <- m+1
        subpop <- colnames(signal[[i]])[1]
        other <- colnames(signal[[i]])[2]
        cell.n <- cellint[cellint$origin==subpop & cellint$destination==other,]
        n.int=c(n.int,paste(subpop, other, sep="-"))
        if (verbose){
          cat("Doing",subpop,"and",other,"...")
        }

        l.both <- intersect(cell.n$ligand[cell.n$origin==subpop],
                          cell.n$ligand[cell.n$origin==other])
        for (l in l.both){
          k <- which(cell.n$ligand==l & cell.n$origin==other)
          cell.n$ligand[k] <- paste("bis",l)
        }
        r.both <- intersect(cell.n$receptor[cell.n$origin==subpop],
                          cell.n$receptor[cell.n$origin==other])
        for (r in r.both){
          k <- which(cell.n$receptor==r & cell.n$origin==other)
          cell.n$receptor[k] <- paste("bis",r)
        }

        if (!is.null(dim(cell.n))){
          g.cell <- graph_from_data_frame(cell.n,directed=TRUE)
          genes <- c(cell.n$ligand,cell.n$receptor)
          prop.genes <- unique(data.frame(gene=genes,display.name=
                                          gsub("^bis ","",genes),
                                        gene.type=c(rep("ligand",nrow(cell.n)),
                                                    rep("receptor",
                                                        nrow(cell.n))),
                                        expressed.in=c(cell.n$origin,
                                                       cell.n$destination),
                                        stringsAsFactors=FALSE))
          prop.genes <- subset(prop.genes,!duplicated(prop.genes[[1]]))
          rownames(prop.genes) <- prop.genes[[1]]
          g.cell <- g.cell %>% set_vertex_attr(
            name="display.name",value=prop.genes[V(g.cell)$name,
                                               "display.name"]) %>%
            set_vertex_attr(name="gene.type",value=prop.genes[V(g.cell)$name,
                                                            "gene.type"]) %>%
            set_vertex_attr(name="expressed.in",value=prop.genes[V(g.cell)$name,
                                                               "expressed.in"])

          if (write){
            write.graph(g.cell,file=paste0('./networks/intercell_network_',
                                         subpop,'-',other,'.graphml'),
                      format="graphml")
          }
          interface[[m]]=cell.n
        } else {
          interface[[m]]="No network"
        }
        cat(' OK', fill=TRUE)
      }
      
    }
    names(interface) <- n.int


    # Full network ----------------
    g.cell <- graph_from_data_frame(cellint[,c(1,2,7,8)],directed=TRUE)
    genes <- c(cellint$ligand,cellint$receptor)
    prop.genes <- unique(data.frame(gene=genes,display.name=
                                      c(cellint$ligand.name,
                                        cellint$receptor.name),
                                    gene.type=c(rep("ligand",nrow(cellint)),
                                                rep("receptor",nrow(cellint))),
                                    stringsAsFactors=FALSE))
    prop.genes <- subset(prop.genes,!duplicated(prop.genes[[1]]))
    rownames(prop.genes) <- prop.genes[[1]]
    g.cell <- g.cell %>% set_vertex_attr(name="display.name",
                                         value=prop.genes[V(g.cell)$name,
                                                          "display.name"]) %>%
      set_vertex_attr(name="gene.type",
                      value=prop.genes[V(g.cell)$name,"gene.type"])
    prop.genes <- unique(data.frame(gene=genes,
                                    display.name=c(cellint$ligand.name,
                                                   cellint$receptor.name),
                                    expressed.in=c(cellint$origin,
                                                   cellint$destination),
                                    stringsAsFactors=FALSE))
    prop.genes <- subset(prop.genes,!duplicated(prop.genes[[1]]))
    rownames(prop.genes) <- prop.genes[[1]]
    g.cell <- g.cell %>% set_vertex_attr(name="expressed.in",
                                         value=prop.genes[V(g.cell)$name,
                                                          "expressed.in"])
    if (write){
      write.graph(g.cell,file=
                    paste0('./networks/full-intercellular-network.graphml'),
                  format="graphml")
    }
    if (plot){
      g.plot <- graph_from_data_frame(cellint,directed=FALSE)
      tmp <- do.call(rbind,strsplit(unique(c(cellint$ligand,cellint$receptor)),split=".",fixed=TRUE))
      cr <- rainbow(length(clusterToStudy))
      names(cr) <- clusterToStudy
      V(g.plot)$vertex.label <- do.call(rbind,strsplit(unique(c(cellint$ligand,cellint$receptor)),split=".",fixed=TRUE))[,2]
      V(g.plot)$label.color <- "black"
      V(g.plot)$color <- cr[tmp[,1]]
      V(g.plot)$shape <- c("circle")
      V(g.plot)$shape[unique(c(cellint$ligand,cellint$receptor)) %in% cellint$receptor] <- c("square")
      V(g.plot)$size <- 10
      E(g.plot)$width <- cellint$LRscore*4
      E(g.plot)$color <- "gray30"

      plot(g.plot,vertex.label=V(g.plot)$vertex.label,main="Intercellular communication network")
      legend("bottomleft",legend=clusterToStudy,fill=cr)
      legend("topleft",legend=c("ligand","receptor"),pch=c(1,0))
    }
  }
  res <- list(interface,cellint[,c(1,2,7,8)])
  names(res) <- c("individual-networks","full-network")
  return(res)
}



#' simplify_interactions
#'
#' @param t the network to be simplified
#' @param lr ligand receptor interactions
#' @param autocrine a logical
#'
#' @return t
#' @export
#'
#' @examples
#' t=data.frame(a.gn=c("CEP63","CEP63"),b.gn=c("MZT2A","DYNC1L2"),
#' type=c("in-complex-with","in-complex-with"))
#' simplify_interactions(t)
simplify_interactions <- function(t,lr=NULL,autocrine=FALSE){

  mergeText <- function(a,b){
    if (length(a)==0)
      b
    else
      if (length(b)==0)
        a
    else
      paste(union(strsplit(a,';')[[1]],strsplit(b,';')[[1]]),collapse=';')
  }


  t$detailed.type <- t$type

  # reduce interaction types
  t$type[t$type %in% c('interacts-with','in-complex-with')] <- 'complex'
  t$type[t$type %in% c('chemical-affects','consumption-controlled-by',
                       'controls-expression-of','controls-phosphorylation-of',
                       'controls-production-of','controls-state-change-of',
                       'controls-transport-of',
                       'controls-transport-of-chemical')] <- 'control'
  t$type[t$type %in% c('catalysis-precedes','reacts-with',
                       'used-to-produce')] <- 'reaction'

  # merge duplicated interactions with same reduced type
  key <- paste(t$a.gn,t$b.gn,t$type,sep='|')
  dk <- which(duplicated(key))
  for (i in dk){
    jj <- setdiff(which(t$a.gn==t[i,"a.gn"] & t$b.gn==t[i,"b.gn"] &
                          t$type==t[i,"type"]),i)
    for (j in jj){
      # t[j,"pmid"] <- mergeText(t[j,"pmid"],t[i,"pmid"])
      t[j,"pathway"] <- mergeText(t[j,"pathway"],t[i,"pathway"])
      t[j,"detailed.type"] <- mergeText(t[j,"detailed.type"],t[i,
                                                               "detailed.type"])
    }
  }
  if (length(dk)>0)
    t <- t[-dk,]

  # promote interaction types in case one interaction is still given for
  # multiple "reduced" interaction types
  key <- paste(t$a.gn,t$b.gn,sep='|')
  dk <- which(duplicated(key))
  for (i in dk){
    jj <- setdiff(which(t$a.gn==t[i,"a.gn"] & t$b.gn==t[i,"b.gn"]),i)
    for (j in jj){
      # t[j,"pmid"] <- mergeText(t[j,"pmid"],t[i,"pmid"])
      t[j,"pathway"] <- mergeText(t[j,"pathway"],t[i,"pathway"])
      t[j,"detailed.type"] <- mergeText(t[j,"detailed.type"],t[i,
                                                               "detailed.type"])
      if (t$type[j]!='control'){
        if (t$type[i]=='control' || t$type[i]=='complex')
          t$type[j] <- t$type[i]
      }
    }
  }
  if (length(dk)>0)
    t <- t[-dk,]

  # eliminate reverse interactions
  l.key <- paste(t$a.gn,t$b.gn,sep='|')
  r.key <- paste(t$b.gn,t$a.gn,sep='|')
  rk <- which(r.key%in%l.key)
  to.remove <- NULL
  for (i in rk){
    j <- setdiff(which(t$a.gn==t[i,"b.gn"] & t$b.gn==t[i,"a.gn"]),i)
    if (j < i){
      # t[j,"pmid"] <- mergeText(t[j,"pmid"],t[i,"pmid"])
      t[j,"pathway"] <- mergeText(t[j,"pathway"],t[i,"pathway"])
      t[j,"detailed.type"] <- mergeText(t[j,"detailed.type"],t[i,
                                                               "detailed.type"])
      if (t$type[j]!='control'){
        if (t$type[i]=='control'){
          t$type[j] <- 'control'
          t$a.gn[j] <- t$a.gn[i]
          t$b.gn[j] <- t$b.gn[i]
        }
        else
          if (t$type[i]=='complex')
            t$type[j] <- 'complex'
      }
      to.remove <- c(to.remove,i)
    }
  }
  if (length(to.remove)>0)
    t <- t[-to.remove,]

  # interactions between 2 ligands or 2 receptors are removed unless 'complex'
  bad <- ((t$a.gn%in%lr$ligand & t$b.gn%in%lr$ligand) |
            (t$a.gn%in%lr$receptor & t$b.gn%in%lr$receptor)) & t$type!='complex'
  t <- t[!bad,]

  if (!autocrine && !is.null(lr)){
    # remove LR pairs
    bad <- (t$a.gn%in%lr$ligand & t$b.gn%in%lr$receptor) |
      (t$a.gn%in%lr$receptor & t$b.gn%in%lr$ligand)
    t <- t[!bad,]
  }

  return(t)
}


