#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("minet")

library(tidyverse)
library(readxl)
library(igraph)
library(dplyr)
library(minet)



computeCLR <- function(matrix){
  #applies CLR to infer network from similarityMatrix
  #
  # Args:   
  #  similarityMatrix: a simmetric matrix with the similarity 
  #                    measurements between the objects.
  #
  #  Returns:
  #    weigthed adjacency matrix of the network of interactions
  #    between the objects computed using the CLR algorithm 
  #    Faith et al. 2007 PLoS Biol 5(1)#applies CLR to infer network from similarityMatrix
  #
  #    Requires minet package  Meyer PE,et ak,  (2008)  BMC Bioinformatics, 9.
  
  clr.m <- minet::clr(matrix)
  return(clr.m)
}


computeCLRv2 <- function(similarityMatrix) {
  #applies CLR to infer network from similarityMatrix
  #
  # Args:   
  #  similarityMatrix: a simmetric matrix with the similarity 
  #                    measurements between the objects.
  #
  #  Returns:
  #    weigthed adjacency matrix of the network of interactions
  #    between the objects computed using the CLR algorithm 
  #    Faith et al. 2007 PLoS Biol 5(1)
  
  # Error handling
  if(!isSymmetric(similarityMatrix)) {
    stop("Non symmetric similarity matrices")
  } 
  
  diag(similarityMatrix) <- NA  #Diagonal elements elminated, no autointeractions
  # z-score transfrom by rows
  zc=apply(similarityMatrix, 2, scale,center=TRUE, scale=TRUE) 
  # negative values are set to zero 
  zc[which(zc<0)] <- 0
  
  # the input matrix is simmetric so the transformation by rows
  # is the transposed matrix 
  zr=t(zc)
  
  # produces likelihood estimates
  clr.z=sqrt(zr^2 + zc^2)
  # add object names( if originally present) 
  if(!is.null(rownames(similarityMatrix))) {
    colnames(clr.z) <- rownames(clr.z) <- rownames(similarityMatrix)
  }
  
  diag(clr.z) <- 0 #set diagonal elements to 0 (instead of NA)
  return(clr.z)
  
}


PCLRC <- function(datamatrix, Niter, frac, rank.thr, probThreshold){
  #performs  PCLRC and imposes a threshol to inferr a networkfrom a data matrix
  #
  # Args:
  #   datamatrix: numeric matrix with the objects measurements.
  #              objects in rows, samples in cols.
  #   Niter: Integer value with the number of iterations. Default is 1000
  #   frac: Fraction of the samples to be considered at each iteration.
  #     Default is 0.75 (75%)            
  #   rank.thr: Fraction of the total predicted interactions to be kept at
  #     each iteration. Default is 0.3 (30% highest scoring interations kept).
  #   probThreshold:  probability threshold for an edge to be considered TRUE
  # Returns:
  #   Weighed adjacency matrix  of the network of interactions
  #   between the objects computed using the PCLRC algorithm 
  #   Saccenti et al 2015 J. Proteome Res., 2015, 14 (2)
  
  
  probabilistNetwork <- getPCLRC(datamatrix, Niter=Niter, frac=frac, rank.thr=rank.thr)
  net <-getAdjByThreshold(probabilistNetwork, threshold=probThreshold)
  return(net)
} 

getPCLRC <- function(datamatrix, Niter, frac, rank.thr){
  #performs the PCLRC to infer probabilistic associations from a data matrix
  #
  # Args:
  #   datamatrix: numeric matrix with the objects measurements.
  #              objects in rows, samples in cols.
  #   Niter: Integer value with the number of iterations. Default is 1000
  #   frac: Fraction of the samples to be considered at each iteration.
  #     Default is 0.75 (75%)            
  #   rank.thr: Fraction of the total predicted interactions to be kept at
  #     each iteration. Default is 0.3 (30% highest scoring interations kept).     
  # Returns:
  #   Weighed adjacency matrix  of the network of interactions
  #   between the objects computed using the PCLRC algorithm 
  #   Saccenti et al 2015 J. Proteome Res., 2015, 14 (2)
  
  
  nsamp=round(dim(datamatrix)[1]*frac) #number of samples per iteration
  Nvaliter=0  #number of valid iterations (no NAs generated)
  
  #this table will store the number of times an interaction was selected
  table <- mat.or.vec(dim(datamatrix)[2], dim(datamatrix)[2])
  
  
  for( it in (1:Niter)){
    #randomly select the samples
    samples=sample(dim(datamatrix)[1], nsamp)
    #build similarity matrix based on correlation
    similarityMatrixSubset <- cor(datamatrix[samples,])
    similarityMatrixSubset <- similarityMatrixSubset^2
    #infeer network for this iteration using CLR
    adjSubset=computeCLR(similarityMatrixSubset)
    
    if(!is.na(sum(adjSubset))){ #valid iteration
      #extract highest links
      out=getHighestLinks(adjSubset, rank.thr)
      #collect output
      table <- table+out
      Nvaliter=Nvaliter+1
    }
  }
  return(table/Nvaliter)
} 

getHighestLinks <- function(matrix, rank.thr, verbose=FALSE)    {
  #Obtains the links with the hightes weigth from an weigthed adjacency matrix
  #
  # Args:
  #   weighedAdjacencyMatrix:  a weighted adjacency matrix (symmetric matrix)
  
  #   rank.thre: fraction of interactions to keep. Default0.3 
  #   Returns: binary matrix with only the higher links having non null values
  
  
  
  if(max(matrix, na.rm=TRUE)==0 ){
    th=0.01
    if(verbose) cat("null clr matrix \n")
  } else{                                      #get the threshold
    if(is.null(rank.thr)) {
      th <- min(matrix[matrix>0], na.rm=TRUE) 
    } else{
      th <- quantile((matrix[upper.tri(matrix)]), 1-rank.thr,na.rm=TRUE) 
      if(th==0){
        if(verbose) cat("threshold too low,  min of the (non-null) matrix chosen instead\n")
        th <- min(matrix[matrix>0])
      }
    }
  }
  ##select the ones that are above the threshold...
  net <- matrix
  net[which(net < th)] =0
  net[which(net >= th)] =1
  return(net)
}                    


getAdjByThreshold <- function(matrix, threshold=0, posDef=TRUE){
  # Impose a threshold on a weigthed adjacency matrix to get connectivity matrix
  #
  # Args:
  #  matrix:  a weighted adjacency matrix (symmetric numeric matrix)
  #  threshold: numeric to be imposed
  #             defaults to zero
  #  posDef: if the matrix is positive define. If not the abs value is considered
  #          default true
  
  x=matrix
  if(posDef) {
    if(sum(x<0)>0) stop(" Negative elments in the net. set posDef  to FALSE" )
    x[which(x>threshold)] <- 1
    x[which(x <= threshold)] <- 0
    
  } else {
    x[which(abs(x)>threshold)] <- 1
    x[which(abs(x) <= threshold)] <- 0
  }
  
  return(x)
}



set_tbl<-read_excel("unique_19_each_585_new.xlsx", sheet="585_1")
set_tbl[,6:100][set_tbl[,6:100 ]==0]<- NA
evo_set_tbl <-set_tbl[,6:24] %>%cbind(set_tbl$cluster)
evoWT_set_tbl<-set_tbl[,25:43] %>%cbind(set_tbl$cluster)
unevo_set_tbl<-set_tbl[,63:81] %>%cbind(set_tbl$cluster)
WT_set_tbl<-set_tbl[,82:100] %>%cbind(set_tbl$cluster)
mut_set_tbl<-set_tbl[,44:62] %>%cbind(set_tbl$cluster)


WT_set_tbl <-WT_set_tbl %>% tibble::column_to_rownames(var="set_tbl$cluster")%>%data.matrix()%>%t()
evo_set_tbl <-evo_set_tbl %>% tibble::column_to_rownames(var="set_tbl$cluster")%>%data.matrix()%>%t()
unevo_set_tbl <-unevo_set_tbl %>% tibble::column_to_rownames(var="set_tbl$cluster")%>%data.matrix()%>%t()
evoWT_set_tbl <-evoWT_set_tbl %>% tibble::column_to_rownames(var="set_tbl$cluster")%>%data.matrix()%>%t()
mut_set_tbl <-mut_set_tbl %>% tibble::column_to_rownames(var="set_tbl$cluster")%>%data.matrix()%>%t()

evo_set_tbl=log(evo_set_tbl)
evoWT_set_tbl=log(evoWT_set_tbl)
unevo_set_tbl=log(unevo_set_tbl)
WT_set_tbl=log(WT_set_tbl)
mut_set_tbl=log(mut_set_tbl)

#Parameter definition
## Niter: Integer value with the number of iterations. Default is 1000 
Niter=1000 
##Fraction of the samples to be considered at each iteration.
#Default is 0.75 (75%)    
frac=0.75

##Fraction of the total predicted interactions to be kept at
##each iteration. Default is 0.3 (30% highest scoring interations kept)
rank.thr=0.3

##Probability threshold for an edge to be considered TRUE
#Default 0.95
probThreshold=0.95

WT_PCLRC <-PCLRC(WT_set_tbl, Niter=Niter, frac=frac, rank.thr=rank.thr ,probThreshold=probThreshold)
WT_PCLRC_net <- graph.adjacency(WT_PCLRC, mode="undirected")
WT_PCLRC_edges <- as.data.frame(get.edgelist(WT_PCLRC_net))%>%`colnames<-`(c("M1","M2"))

unevo_PCLRC <-PCLRC(unevo_set_tbl, Niter=Niter, frac=frac, rank.thr=rank.thr ,probThreshold=probThreshold)
unevo_PCLRC_net <- graph.adjacency(unevo_PCLRC, mode="undirected")
unevo_PCLRC_edges <- as.data.frame(get.edgelist(unevo_PCLRC_net))%>%`colnames<-`(c("M1","M2"))

evo_PCLRC <-PCLRC(evo_set_tbl, Niter=Niter, frac=frac, rank.thr=rank.thr ,probThreshold=probThreshold)
evo_PCLRC_net <- graph.adjacency(evo_PCLRC, mode="undirected")
evo_PCLRC_edges <- as.data.frame(get.edgelist(evo_PCLRC_net))%>%`colnames<-`(c("M1","M2"))


evoWT_PCLRC <-PCLRC(evoWT_set_tbl, Niter=Niter, frac=frac, rank.thr=rank.thr ,probThreshold=probThreshold)
evoWT_PCLRC_net <- graph.adjacency(evoWT_PCLRC, mode="undirected")
evoWT_PCLRC_edges <- as.data.frame(get.edgelist(evoWT_PCLRC_net)%>%`colnames<-`(c("M1","M2")))

mut_PCLRC <-PCLRC(mut_set_tbl, Niter=Niter, frac=frac, rank.thr=rank.thr ,probThreshold=probThreshold)
mut_PCLRC_net <- graph.adjacency(mut_PCLRC, mode="undirected")
mut_PCLRC_edges <- as.data.frame(get.edgelist(mut_PCLRC_net))%>%`colnames<-`(c("M1","M2"))


write.csv(WT_PCLRC_edges,"WT_PCLRC_edges_1.csv")
write.csv(unevo_PCLRC_edges,"unevo_PCLRC_edges_1.csv")
write.csv(evo_PCLRC_edges,"evo_PCLRC_edges_1.csv")
write.csv(evoWT_PCLRC_edges,"evoWT_PCLRC_edges_1.csv")
write.csv(mut_PCLRC_edges,"mut_PCLRC_edges_1.csv")



set_tbl<-read_excel("unique_19_each_585_new.xlsx", sheet="585_2")
set_tbl[,6:100][set_tbl[,6:100 ]==0]<- NA
evo_set_tbl <-set_tbl[,6:24] %>%cbind(set_tbl$cluster)
evoWT_set_tbl<-set_tbl[,25:43] %>%cbind(set_tbl$cluster)
unevo_set_tbl<-set_tbl[,63:81] %>%cbind(set_tbl$cluster)
WT_set_tbl<-set_tbl[,82:100] %>%cbind(set_tbl$cluster)
mut_set_tbl<-set_tbl[,44:62] %>%cbind(set_tbl$cluster)


WT_set_tbl <-WT_set_tbl %>% tibble::column_to_rownames(var="set_tbl$cluster")%>%data.matrix()%>%t()
evo_set_tbl <-evo_set_tbl %>% tibble::column_to_rownames(var="set_tbl$cluster")%>%data.matrix()%>%t()
unevo_set_tbl <-unevo_set_tbl %>% tibble::column_to_rownames(var="set_tbl$cluster")%>%data.matrix()%>%t()
evoWT_set_tbl <-evoWT_set_tbl %>% tibble::column_to_rownames(var="set_tbl$cluster")%>%data.matrix()%>%t()
mut_set_tbl <-mut_set_tbl %>% tibble::column_to_rownames(var="set_tbl$cluster")%>%data.matrix()%>%t()

evo_set_tbl=log(evo_set_tbl)
evoWT_set_tbl=log(evoWT_set_tbl)
unevo_set_tbl=log(unevo_set_tbl)
WT_set_tbl=log(WT_set_tbl)
mut_set_tbl=log(mut_set_tbl)



WT_PCLRC <-PCLRC(WT_set_tbl, Niter=Niter, frac=frac, rank.thr=rank.thr ,probThreshold=probThreshold)
WT_PCLRC_net <- graph.adjacency(WT_PCLRC, mode="undirected")
WT_PCLRC_edges <- as.data.frame(get.edgelist(WT_PCLRC_net))%>%`colnames<-`(c("M1","M2"))

unevo_PCLRC <-PCLRC(unevo_set_tbl, Niter=Niter, frac=frac, rank.thr=rank.thr ,probThreshold=probThreshold)
unevo_PCLRC_net <- graph.adjacency(unevo_PCLRC, mode="undirected")
unevo_PCLRC_edges <- as.data.frame(get.edgelist(unevo_PCLRC_net))%>%`colnames<-`(c("M1","M2"))

evo_PCLRC <-PCLRC(evo_set_tbl, Niter=Niter, frac=frac, rank.thr=rank.thr ,probThreshold=probThreshold)
evo_PCLRC_net <- graph.adjacency(evo_PCLRC, mode="undirected")
evo_PCLRC_edges <- as.data.frame(get.edgelist(evo_PCLRC_net))%>%`colnames<-`(c("M1","M2"))


evoWT_PCLRC <-PCLRC(evoWT_set_tbl, Niter=Niter, frac=frac, rank.thr=rank.thr ,probThreshold=probThreshold)
evoWT_PCLRC_net <- graph.adjacency(evoWT_PCLRC, mode="undirected")
evoWT_PCLRC_edges <- as.data.frame(get.edgelist(evoWT_PCLRC_net)%>%`colnames<-`(c("M1","M2")))

mut_PCLRC <-PCLRC(mut_set_tbl, Niter=Niter, frac=frac, rank.thr=rank.thr ,probThreshold=probThreshold)
mut_PCLRC_net <- graph.adjacency(mut_PCLRC, mode="undirected")
mut_PCLRC_edges <- as.data.frame(get.edgelist(mut_PCLRC_net)%>%`colnames<-`(c("M1","M2")))


write.csv(WT_PCLRC_edges,"WT_PCLRC_edges_2.csv")
write.csv(unevo_PCLRC_edges,"unevo_PCLRC_edges_2.csv")
write.csv(evo_PCLRC_edges,"evo_PCLRC_edges_2.csv")
write.csv(evoWT_PCLRC_edges,"evoWT_PCLRC_edges_2.csv")
write.csv(mut_PCLRC_edges,"mut_PCLRC_edges_2.csv")


