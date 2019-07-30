#' scCLUE: single cell clustering through ensemble feature selection and similarity measurement
#' @param log_scdata (M x N) dimensional matrix for gene expressions (log scale)
#' M: The number of genes
#' N: The number of cells
#' @param nCls The number of clusters
#' @param nEns The number of trials for ensemble feature selection and similarity measurements
#' @param K The number of nearest neighbors for KNN networks
#' @param smin Ratio of the minimum feature sampling
#' @param smax Ratio of the maximum feature sampling
#'
#' @author Hyundoo Jeong
#'
#' @references
#' Hyundoo Jeong
#'
#' Effective single cell clustering through ensemble feature selection and similarity measurements
#'
#' @example
#' library(scCLUE)
#' library(edgeR)
#' data(scdata)
#'
#' logsc <- log10(1+cpm(scdata))
#' final_cls <- scCLUE(log_scdata = logsc, nCls = length(table(label)), nEns = 15, K = c(5, 10), nPCs=15)

scCLUE <- function(log_scdata = data, nCls, nEns = 15, K = c(5,10), nPCs=15, smin = 0.6, smax = 1){


  my_features <- IdentifyFeatures(inData = log_scdata, pFeat = 0.2)


  nCells <- dim(log_scdata)[2]
  ensNet <- Matrix(0, nrow = nCells, ncol = nCells, sparse = TRUE)


  for(ix in 1:nEns){
    for(ik in K){
      sample_rng <- round(runif(n=1, min = smin, max = smax), digit = 2)
      sampled_genes <- sample(my_features$VarGenes, size = ceiling(sample_rng*length(my_features$VarGenes)), replace = F)

      umap_res = umap(t(log_scdata[sampled_genes,]), method ='umap-learn')
      umap_knn <- MakeKNN(data = umap_res$layout, k_num = ik)
      umap_knn <- umap_knn + t(umap_knn)


      my_pca_log <- prcomp_irlba(t(log_scdata[sampled_genes,]), n = nPCs)
      PCs <- (my_pca_log$x)
      cor_dist <- cor(t(PCs), method ='pearson')
      cor_knn <- make.kNNG(1-cor_dist, k = ik, symm = F)
      cor_knn <- cor_knn + t(cor_knn)

      temp <- umap_knn + cor_knn
      ensNet <- ensNet + make_dsm(temp)
    }
  }


  gg <- graph_from_adjacency_matrix(ensNet, mode='undirected', weighted =T)
  cl <- cluster_louvain(gg)


  sind <- which(table(cl$membership) == 1)
  if(length(sind) != 0){
    clid <- setdiff(names(table(cl$membership)), names(sind))
    for(ss in sind){
      single <- which(cl$membership == ss)
      mscore <- matrix(0, ncol = length(clid), nrow=1)
      ii <- 1
      for(cc in clid){
        cid <- which(cl$membership == cc)
        mscore[ii] <- mean(ensNet[single, cid])
        ii <- ii+1
      }
      cl$membership[single] <- clid[which(mscore == max(mscore))]
    }

    clid <- names(table(cl$membership))
    new_mem <- cl$membership
    cid <- 1
    for(mk in 1:length(clid)){
      idx <- which(cl$membership == clid[mk])
      new_mem[idx] <- mk
      cid <- cid +1
    }
    cl$membership <- new_mem
  }


  sep_scores <- ComputeSepScores(log_scdata = log_scdata, membership= cl$membership, topK = my_features$topK)

  final_cls<- MergeClusters(cl$membership, sep_scores)
  id <- which(names(final_cls) == nCls)
  if(length(id) == 0){
    id <- max(names(final_cls))
  }
  final_cls$membership <- final_cls[[id]]


  return(final_cls)
}






#' Identify feature genes (highly variable genes)
#' @param inData The gene expression matrix in log scale. Each row corresponds to the gene and column corresponds to the cell
#' @param pFeat Percentage of potential features to be selected
#'
IdentifyFeatures <- function(inData = data, pFeat = 0.2){


  rm <- rowMeans(inData)
  rv <- rowVars(inData)

  # first filtering based on the mean expression
  hind <- which(rm >= median(rm))
  my_rv <- rv[hind]
  names(my_rv) <- names(hind)
  my_vargenes <- names(head(x = sort(my_rv, decreasing = T), ceiling(pFeat*dim(inData)[1])) )

  topK <- names(head(x = sort(my_rv, decreasing = T), ceiling(0.05*dim(inData)[1])) )


  features <- list("VarGenes" = my_vargenes,
                   "topK" = topK

  )
  return(features)
}





#' Construct a KNN (K nearest neighbors) network
#' @param data The gene expression matrix in log scale. Each row corresponds to the gene and column corresponds to the cell
#' @param k_num The number of nearest neighbors to construct KNN network
#'
MakeKNN <- function(data = in_data, k_num = in_k){
  edist <- parDist(data, method = 'euclidean')
  edist <- as.matrix(edist)

  my_knn <- make.kNNG(edist, k = k_num, symm = F)
  my_knn <- Matrix(my_knn, sparse = TRUE)
  return(my_knn)
}






#' Compute a separation score between different clusters
#'
#' ComputeSepScores computes the separation score between different clsters using Welch's t-test
#'
#' @param data The gene expression matrix in log scale. Each row corresponds to the gene and column corresponds to the cell
#' @param membership Membership index for the clustering results
#' @param topK Top K highly variable genes
#'
ComputeSepScores <- function(log_scdata = in_data, membership, topK){

  sep_scores <- matrix(Inf, nrow = length(table(membership)),  ncol = length(table(membership)))
  for(ix in 1:(length(table(membership))-1)){
    for(iy in (ix+1):length(table(membership))){

      pva <- matrix( nrow = length(topK), ncol=1)
      id1 <- which(membership == ix)
      id2 <- which(membership == iy)

      for(i in 1:length(topK)){
        a <- t.test(log_scdata[topK[i], id1] , log_scdata[topK[i], id2])
        pva[i] <- a$p.value
      }

      nan.id <- which(pva == 'NaN')
      pva[nan.id] <- Inf

      sep_scores[ix, iy] <- -log10(min(pva))

    }
  }
  return(sep_scores)
}






#' Merge different clusters based on a separation score
#'
#' MergeClusters iteratively merges different clusters based on the separation scores.
#'
#' @param membership Membership index for the clustering results
#' @param sep_scores Separation scores between each clusters
MergeClusters <- function(membership, sep_scores){

  ncl <- length(table(membership))
  my_membership <- list()
  my_membership[[toString(ncl)]] <- membership
  cl_id <- as.numeric(names(table(membership)))

  mid <- 1
  cnum <- ncl
  for( ic in (ncl-1):1){
    cur_mem <- my_membership[[mid]]
    mem_id <- as.numeric(names(table(cur_mem)))

    ind <- which(sep_scores == min(sep_scores[mem_id, mem_id]), arr.ind = TRUE)
    my_id <- which(cur_mem == cl_id[max(ind)])


    cur_mem[my_id] <- cl_id[min(ind)]

    mid <- mid + 1
    my_membership[[toString(ic)]] <- cur_mem
  }


  new_membership <- list()
  for(ic in 1:length(my_membership)){
    memid <- names(table(my_membership[[ic]]))

    new_mem <- my_membership[[ic]]
    for (k in 1:length(memid)){
      ids <- which(my_membership[[ic]] == memid[k])
      new_mem[ids] <- k
    }
    new_membership[[toString(length(memid))]] <- new_mem
  }


  return(new_membership)


}



#' Make a doubly stochastic matrix
#'
#' make_dsm converts a square matrix into a doubly stochastic matrix.
#'
#' @param my_may Square matrix
make_dsm <- function(my_mat){
  D1 <- Matrix(diag(1/Matrix::rowSums(my_mat)), sparse = TRUE)
  TT <- D1 %*% my_mat
  D2 <- Matrix(diag(1/sqrt(Matrix::colSums(TT))), sparse = TRUE)
  TT <- TT %*% D2
  my_dsm <- TT %*% t(TT)
  return(my_dsm)

}


