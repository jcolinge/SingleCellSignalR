#' @title intra network
#' @description Computes intracellular networks linked to genes of interest.
#'
#' @details
#' coi must correspond to one of the object's cluster names. If no names was given
#' the clusters are numbered such as "cluster 1", "cluster 2", ...
#' @details
#' `cell.prop` is set to 0.2
#' by default to avoid unreadable downstream networks. However if the calculated
#' network is too small or non-existent (or too big) the user can try lower
#' (or higher) values.
#' @details
#' In the case of mouse data, the function converts mouse genes
#' in human orthologs (according to Ensembl) such that the Reactome/KEGG
#' interaction database can be exploited, and finally output genes are converted
#' back to mouse.
#' @details
#' If `write` is TRUE, then the function writes two different
#' files. A graphML file in the *network* folder for intracellular
#' interactions downstream the gene of interest (goi) named
#' "intracell_network_coi-receptors.graphml". A text file in the *network*
#' folder containing the information about the pathways in which the interactions
#' are in, named "intracell_network_pathway_analysis_coi-receptors.txt".
#'
#' @param dm a SCSRDataModel object
#' @param goi gene of interest (typically a receptor)
#' @param coi name of the cluster of interest
#' @param obj a SCSRInference object
#' @param cell.prop a threshold, only the genes expressed in this proportion of
#' the cells of the coi will be taken into account
#' @param write a logical (if TRUE writes graphML and text files for the
#' internal networks)
#' @param plot a logical
#' @param add.lig a logical (if TRUE adds the goi associated ligands from
#' signal to the network)
#' @param connected a logical (if TRUE keeps only the genes connected to
#' the goi)
#' @param verbose a logical
#'
#' @return The function returns a list containing the internal networks
#' linked to the genes of interest (goi)
#'
#' @export
#'
#' @importFrom foreach foreach %do%
#' @importFrom stats phyper
#' @importFrom multtest mt.rawp2adjp
#' @importFrom utils write.table
#' @import igraph
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
#' net <- intra_network(obj, "gene 20", 1, obj.int)

intra_network <- function(dm,goi,coi,obj=NULL,cell.prop=0.2,
                         write=TRUE,plot=TRUE,add.lig=TRUE,
                         most.variables=TRUE,connected=FALSE,verbose=TRUE){
  if (!is(dm, "SCSRDataModel")){
        stop("dm must be a SCSRDataModel object")
    }
    if (!is.null(obj)){
    signal <- obj@LRinter
      if (!is(obj, "SCSRInference")){
        stop("obj must be a SCSRInference object")
      }
    }

  c.names <- dm@cluster$names
  cluster <- dm@cluster$id
  data <- dm@ncounts$matrix
  genes <- rownames(dm@ncounts$matrix)

  if (!is.null(dm@ncounts$matrix.mv)&most.variables){
        cat("Matrix of most variable genes used. To use the whole matrix set most.variables 
            parameter to FALSE.\n")
        data <- dm@ncounts$matrix.mv
        genes <- rownames(dm@ncounts$matrix.mv)
    }
  if (dir.exists("networks")==FALSE & write==TRUE){
    dir.create("networks")
  }
  if (is.null(c.names)==TRUE){
    c.names <- paste("cluster",seq_len(max(cluster)))
  }
  if (min(cluster)!=1){
    cluster <- cluster + 1 - min(cluster)
  }
  if (length(c.names)!=max(cluster) | sum(duplicated(c.names))>0 |
      grepl("/",paste(c.names,collapse =""))){
    stop("The length of c.names must be equal to the number of clusters and must
        contain no duplicates. The cluster names must not include special
        characters")
  }
  if (!is.element(coi,c.names)){
    stop(paste(coi,"must be included in c.names.","If c.names is not provided, 
      it is set to cluster 1, cluster 2, ...,
        cluster N. WIth N the maximum number of clusters"))
  }
  opar <- par()
  species <- dm@initial.organism
  goi.ini <- goi
  if (species!='hsapiens'){
    ortho <- data.frame(Hsapiens = rownames(data),
              species = object@ncounts$initial.orthologs)
    goi.ini <- goi
    goi <- as.character(ortho[ortho$species == goi, "Hsapiens"])
  }
  pw.names <- strsplit(PwC_ReactomeKEGG$pathway,";")
  pw.sizes <- table(unlist(pw.names))
  max.pw.size <- 500
  good.pw <- pw.sizes[pw.sizes<=max.pw.size]

  data.tmp <- data[,cluster==which(c.names==coi)]
  data.tmp <- data.tmp[rowSums(data.tmp)>0,]

  # expressed genes
  good <- apply(data.tmp,1,function(x) sum(x>0)/ncol(data.tmp)>cell.prop)
  visible.genes <- unique(c(rownames(data.tmp)[good],goi))
  visible.n <- PwC_ReactomeKEGG[PwC_ReactomeKEGG$a.gn%in%visible.genes &
                                 PwC_ReactomeKEGG$b.gn%in%visible.genes,]
  red.visible.n <- simplify_interactions(visible.n,LRdb)
  res <- list()
  qq <- 0

  for (receptors in goi){
    qq <- qq+1
    if (!is.element(receptors,rownames(data.tmp))){
      cat(paste0(goi.ini[qq]," is not expressed in ",coi),fill=TRUE)
    } else {
      # receptor containing pathways ------------
      contains.receptors <- intersect(unlist(
        pw.names[PwC_ReactomeKEGG$a.gn%in%receptors |
                   PwC_ReactomeKEGG$b.gn%in%receptors]),names(good.pw))
      if (verbose){
        cat(paste0("Patwhay(s) that include ",goi.ini[qq],":"),fill=TRUE)
        for (i in contains.receptors){
          cat(paste0("   - ",i),fill=TRUE)
        }
      }
      if (sum(grepl("added",contains.receptors))>1){
        contains.receptors <- contains.receptors[contains.receptors!="added"]
      }
      contain.n=NULL
      for (i in contains.receptors){
        contain.n <- rbind(contain.n,PwC_ReactomeKEGG[
          grepl(i,PwC_ReactomeKEGG$pathway),])
      }
      contain.n <- unique(rbind(contain.n,PwC_ReactomeKEGG[
        PwC_ReactomeKEGG$a.gn%in%visible.genes &
          PwC_ReactomeKEGG$b.gn%in%receptors,],
        PwC_ReactomeKEGG[PwC_ReactomeKEGG$a.gn%in%receptors &
                           PwC_ReactomeKEGG$b.gn%in%visible.genes,]))
      red.contain.n <- simplify_interactions(contain.n)

      # intersect corr network and receptor-containing network ------------
      key.visible <- paste(red.visible.n$a.gn,red.visible.n$b.gn,sep="|")
      key.contain <- paste(red.contain.n$a.gn,red.contain.n$b.gn,sep="|")
      net.n <- red.visible.n[key.visible%in%key.contain,]

      if (nrow(net.n)>0){
        # add ligands
        add.net <- NULL
        nam=NULL
        siz=NULL
        if (!is.null(signal) & add.lig){
          for (i in names(signal)){
            if (grepl(paste0(coi,"$"),i)){
              tmp=signal[[i]]
              if (species!="hsapiens"){
                ortho <- data.frame(Hsapiens = rownames(data),
                                        species = object@ncounts$initial.orthologs)
                lig <- data.frame(species = tmp$ligand, id = seq(1, nrow(tmp), 1))
                rec <- data.frame(species = tmp$receptor, id = seq(1, nrow(tmp), 1))
                lig <- merge(lig, ortho, by.x = "species", order = FALSE)
                lig <- lig[order(lig$id),]
                rec <- merge(rec, ortho, by.x = "species", order = FALSE)
                rec <- rec[order(rec$id),]
                tmp$ligand <- lig$Hsapiens
                tmp$receptor <- rec$Hsapiens
              }
              if (is.element(receptors,tmp[,2])){
                if(nrow(data[tmp[tmp[,2]%in%receptors,1],
                                          cluster==as.numeric(which(
                                            c.names%in%colnames(tmp)[1]))])>1){
                  siz <- c(siz,rowMeans(data[tmp[tmp[,2]%in%receptors,1],
                                          cluster==as.numeric(which(
                                            c.names%in%colnames(tmp)[1]))]))
                  }else{
                    siz <- c(siz,data[tmp[tmp[,2]%in%receptors,1],
                                          cluster==as.numeric(which(
                                            c.names%in%colnames(tmp)[1]))])
                  }
                
                nam.tmp <- rep(colnames(tmp)[1],nrow(tmp[tmp[,2]%in%receptors,]))
                colnames(tmp)[seq_len(2)] <- c("a.gn","b.gn")
                tmp[,1] <- paste(nam.tmp[1],tmp[,1],sep="-")
                add.net <- rbind(add.net,tmp[tmp[,2]%in%receptors,])
                nam <- c(nam,nam.tmp)
              }
            }
          }
          if (is.null(add.net)==FALSE){
            add.net <- cbind(add.net[,seq_len(2)],nam,location="extra",
                            type="control")
          }
        }
        net.tmp <- cbind(net.n[,seq_len(2)],nam=rep(coi,nrow(net.n)),
                        location="intra",type=net.n$type)
        net.f <- rbind(add.net,net.tmp)
        if (species!="hsapiens"){
          lig <- data.frame(Hsapiens = net.f[,1], id = seq(1, nrow(net.f), 1))
          rec <- data.frame(Hsapiens = net.f[,2], id = seq(1, nrow(net.f), 1))
          lig <- merge(lig, ortho, by.x = "Hsapiens", order = FALSE)
          lig <- lig[order(lig$id),]
          rec <- merge(rec, ortho, by.x = "Hsapiens", order = FALSE)
          rec <- rec[order(rec$id),]
          net.f[,1] <- lig$species
          net.f[,2] <- rec$species
        }
        g.net <- graph_from_data_frame(net.f,directed=TRUE)
        g.net.tmp <- as.undirected(g.net)

        y <- shortest_paths(g.net.tmp,receptors,V(g.net.tmp))
        g.net.tmp <- graph_from_data_frame(net.f[,seq_len(2)],directed=FALSE)
        V(g.net.tmp)$status <- "pw.related"
        V(g.net.tmp)$status[
          which(unique(c(net.f$a.gn,net.f$b.gn)) %in% receptors)] =
          "gene.of.interest"
        V(g.net.tmp)$status[which(unique(c(net.f$a.gn,net.f$b.gn)) %in%
                                    net.f$a.gn[net.f$location=="extra"])] =
          "ligand"
        E(g.net.tmp)$int.type <- as.character(net.f$type)

        y <- unlist(lapply(y$vpath,function(x) length(x)))
        y <- y-1
        if (sum(y==-1)>0 & !connected){
          y[y==-1] <- 1
        }
        if (sum(y==-1)>0 & connected){
          nam.tmp <- unique(c(net.f$a.gn,net.f$b.gn))[y==-1]
          net.f <- net.f[!net.f$a.gn%in%nam.tmp & !net.f$b.gn%in%nam.tmp,]
          y <- y[y!=-1]
          g.net <- graph_from_data_frame(net.f,directed=TRUE)
        }
        names(y) <- unique(c(net.f$a.gn,net.f$b.gn))
        y[net.f$a.gn[net.f$location=="extra"]]=-1
        x=vector("numeric",length=length(y))
        names(x) <- unique(c(net.f$a.gn,net.f$b.gn))
        if (sum(net.f$location=="extra")){
          x[net.f$a.gn[net.f$location=="extra"]] <- (seq(0,length(
            net.f$a.gn[net.f$location=="extra"])*4-1,4)-max(seq(0,length(
              net.f$a.gn[net.f$location=="extra"])*4-1,4))/2)
        }
        x[receptors] <- 0
        for (j in seq_len(max(y))){
          if (sum(y==j)>1){
            x[y==j] <- seq(-10,10,19/sum(y==j))[seq_len(sum(y==j))]
          }
        }
        l <- cbind(x,-y)
        for (i in seq_len(max(y))){
          if(sum(y==i)>1){
            index <- which(y==i)
            l[index[seq(1,length(index),2)],2] <-
              l[index[seq(1,length(index),2)],2]-0.2
            l[index[seq(2,length(index),2)],2] <-
              l[index[seq(2,length(index),2)],2]+0.2
          }
        }

        rownames(l) <- unique(c(net.f$a.gn,net.f$b.gn))

        if (sum(net.f$location %in% "intra")==0){
          V(g.net)$size[y==-1] <- (log(siz)+5)*4
          V(g.net)$size[V(g.net)$size<exp(-5)] <- exp(-5)
        } else {
          V(g.net)$size <- (log(c(siz,rowMeans(data.tmp[
            unique(c(net.n$a.gn,net.n$b.gn)),])))+5)*4
          V(g.net)$size[V(g.net)$size<exp(-5)]=exp(-5)
        }
        V(g.net)$size[which(unique(c(net.f$a.gn,net.f$b.gn)) %in% receptors)] <-
          35
        V(g.net)$shape <- c("circle")
        V(g.net)$shape[which(unique(c(net.f$a.gn,net.f$b.gn)) %in% receptors)] <-
          c("rectangle")
        V(g.net)$color <- c("lightcyan")
        V(g.net)$color[which(unique(c(net.f$a.gn,net.f$b.gn)) %in% receptors)] <-
          c("indianred1")
        V(g.net)$color[which(unique(c(net.f$a.gn,net.f$b.gn)) %in%
                               net.f$a.gn[net.f$location=="extra"])] <-
          c("lightyellow1")
        V(g.net)$label.color <- "black"
        V(g.net)$label.dist[which(!unique(c(net.f$a.gn,net.f$b.gn))
                                  %in% receptors)] <- 1
        for (i in unique(l[,2])){
          if (i!=0){
            V(g.net)$label.degree[l[,2]==i] <- rep(c(pi/2,-pi/2),length(y))[
              seq_len(sum(l[,2]==i))]
          }
        }
        V(g.net)$label.dist[which(unique(c(net.f$a.gn,net.f$b.gn)) %in%
                                    receptors)] <- 0
        V(g.net)$label.degree[which(unique(c(net.f$a.gn,net.f$b.gn)) %in%
                                      receptors)] <- 0
        E(g.net)$arrow.mode <- as.numeric(grepl("control",net.f$type))*2
        E(g.net)$arrow.size <- rep(0.4,nrow(net.f))
        E(g.net)$color <- "gray"

        # pathway enrichment statistics
        g.pw.names <- strsplit(net.n$pathway,";")
        g.pw.names <- g.pw.names[g.pw.names!="added"]
        g.pw.sizes <- table(unlist(g.pw.names))
        N <- nrow(PwC_ReactomeKEGG)
        n <- nrow(net.n)
        pw.table <- foreach(pw=names(g.pw.sizes),.combine=rbind) %do% {
          K <- pw.sizes[pw]
          if (K <= max.pw.size){
            k <- g.pw.sizes[pw]
            pval <- 1-phyper(q=k-1,m=K,n=N-K,k=n)
            data.frame(pathway=pw,in.pw=k,pw.size=K,pval=pval,
                       stringsAsFactors=FALSE)
          }
          else
            NULL
        }

        par(mfrow=c(2,1))
        par(las=2)
        par(mar=c(0,0,1,0))

        if (species!="hsapiens"){
          lab = unique(c(net.f$a.gn,net.f$b.gn))

          li = lab[net.f$location=="extra"]
          li = do.call(rbind,strsplit(li,split = "-"))

          change <- data.frame(Hsapiens = li[,2], id = seq(1, nrow(li), 1))
          change <- merge(change, ortho, by.x = "Hsapiens", order = FALSE)
          change <- change[order(change$id),]
          li[,2] <- change$species
          
          li = paste(li[,1],li[,2],sep="-")

          ot = lab[net.f$location!="extra"]
          change <- data.frame(Hsapiens = ot, id = seq(1, length(ot), 1))
          change <- merge(change, ortho, by.x = "Hsapiens", order = FALSE)
          change <- change[order(change$id),]
          ot <- change$species

          plot(g.net, layout=l,main = coi,vertex.label = c(li,ot))

        } else {
          plot(g.net, layout=l,main = coi)
        }

        if (!is.null(pw.table)){
          rawp <- pw.table$pval
          if (length(rawp)==1){
            adj <- rawp
            qval <- adj
          } else {
            adj <- mt.rawp2adjp(rawp,"BH")
            qval <- adj$adjp[order(adj$index),"BH"]
          }

          pw.table <- cbind(pw.table,data.frame(qval=qval))
          rownames(pw.table) <- NULL
          pw.table <- pw.table[pw.table$qval<0.05,]
          if (nrow(pw.table)>0){
            pw.table <- pw.table[order(pw.table$in.pw,decreasing=FALSE),]
            nc <- max(nchar(as.character(pw.table$pathway)))*0.3

            mtext(paste(receptors,"related pathways"),side=1,line=2,
                  las=FALSE)
            par(mar=c(2,nc,4,2))
            barplot((pw.table$in.pw),horiz=TRUE,
                    names.arg=(paste(pw.table$pathway,"*")),
                    cex.names=0.7,col="azure2",border="gray60")
          } else {
            if (verbose){
              cat("No significant associated pathway",fill=TRUE)
            }
          }
        } else {
          if (verbose){
            cat(paste("No associated genes downstream",receptors, "in",coi),
                fill=TRUE)
          }
        }

        if (write){
          if (!is.null(pw.table)){
            write.table(pw.table[order(pw.table$qval),],
                        file=paste0(
                          './networks/intracell_network_pathway_analysis_',
                                    coi,'-',receptors,'.txt'),sep="\t",
                        quote=FALSE,row.names=FALSE)
          }
          write.graph(g.net.tmp,paste0('./networks/intracell_network_',coi,
                                       '-',receptors,'.graphml'),
                      format="graphml")
        }
        res[[qq]] <- net.f
      }
      else {
        cat("No interactions found.", fill=TRUE)
      }
    }
  }
  par(opar)
  return(res)
}
