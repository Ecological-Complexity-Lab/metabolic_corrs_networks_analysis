
.libPaths("/gpfs0/shai/projects/R4/R-4.2.0/lib64/R/library")
print(.libPaths())
print(sessionInfo())
Sys.setlocale(category = "LC_ALL", locale = "")


library(tidyverse)
library(ggplot2)
library(igraph)
library(data.table)
library(dplyr)

common_WT_pos <- read.csv("common_WT_pos_PCLRC.csv",row.names = 1)
common_WT_neg <- read.csv("common_WT_neg_PCLRC.csv",row.names = 1)

common_evoWT_pos <- read.csv("common_evoWT_pos_PCLRC.csv",row.names = 1)
common_evoWT_neg <- read.csv("common_evoWT_neg_PCLRC.csv",row.names = 1)

common_evo_pos <- read.csv("common_evo_pos_PCLRC.csv",row.names = 1)
common_evo_neg <- read.csv("common_evo_neg_PCLRC.csv",row.names = 1)

common_unevo_pos <- read.csv("common_unevo_pos_PCLRC.csv",row.names = 1)
common_unevo_neg <- read.csv("common_unevo_neg_PCLRC.csv",row.names = 1)

common_mut_pos <- read.csv("common_mut_pos_PCLRC.csv",row.names = 1)
common_mut_neg <- read.csv("common_mut_neg_PCLRC.csv",row.names = 1)

print("finished_common")

create_igraph <- function(df){
  before_igraph <- df%>% dplyr::select(from=M1, to=M2, weight=cor)
  igraph <-igraph::graph.data.frame(before_igraph ,directed=F) }

evo_negative_p_igraph_df <- create_igraph(common_evo_neg)
evo_positive_p_igraph_df <-create_igraph(common_evo_pos) 
evo=list(evo_positive_p_igraph_df,evo_negative_p_igraph_df)

unevo_negative_p_igraph_df <-create_igraph(common_unevo_neg) 
unevo_positive_p_igraph_df <-create_igraph(common_unevo_pos) 
unevo=list(unevo_positive_p_igraph_df,unevo_negative_p_igraph_df)

WT_negative_p_igraph_df <- create_igraph(common_WT_neg)
WT_positive_p_igraph_df <- create_igraph(common_WT_pos)
WT=list(WT_positive_p_igraph_df,WT_negative_p_igraph_df)

evoWT_negative_p_igraph_df <- create_igraph(common_evoWT_neg)
evoWT_positive_p_igraph_df <-create_igraph( common_evoWT_pos) 
evoWT=list(evoWT_positive_p_igraph_df,evoWT_negative_p_igraph_df)

mut_negative_p_igraph_df <- create_igraph(common_mut_neg)
mut_positive_p_igraph_df <-create_igraph(common_mut_pos) 
mut=list(mut_positive_p_igraph_df,mut_negative_p_igraph_df)

print("finished_igraph")

#uploading the random networks

WT_shuffle_networks <- readRDS("shuffled_nets/WT_shuffle_networks_pos_neg_PCLRC.RData")[1:500]
evo_shuffle_networks <- readRDS("shuffled_nets/evo_shuffle_networks_pos_neg_PCLRC.RData")[1:500]
unevo_shuffle_networks <- readRDS("shuffled_nets/unevo_shuffle_networks_pos_neg_PCLRC.RData")[1:500]
evoWT_shuffle_networks <- readRDS("shuffled_nets/evoWT_shuffle_networks_pos_neg_PCLRC.RData")[1:500]
mut_shuffle_networks <- readRDS("shuffled_nets/mut_shuffle_networks_pos_neg_PCLRC.RData")[1:500]

print("finished_Rdata")

jaccard_rnd_vs_real_networks <- function(list_of_rnd, list_of_true){
  res <- c()
  bact <- c()
  for (set_2_num in c(1:length(list_of_true))){ #1:5- true
    for (set_1_num in c(1:length(list_of_rnd))){ #1:5- rnd
      if (set_2_num!=set_1_num){
        ran_res <- c()
        for (ran_net in c(1:length(list_of_rnd[[set_1_num]]))){ #1:500
          rnd=list_of_rnd[[set_1_num]][[ran_net]]
          g_pos <- delete.edges(rnd, which(E(rnd)$weight <0))
          g_neg <- delete.edges(rnd, which(E(rnd)$weight >0))
          pos_common=length(E(g_pos %s% list_of_true[[set_2_num]][[1]]))
          neg_common=length(E(g_neg %s% list_of_true[[set_2_num]][[2]]))
          common= pos_common+neg_common
          
          bact <- append(bact,paste0("rnd_",set_1_num,"_","true_",set_2_num))
          print(paste0("pos_common=",pos_common,"_neg_common=",neg_common))
          print(paste0("rnd_pos=",length(E(g_pos)),"_rnd_neg=",length(E(g_neg)),"_true_pos=",length(E(list_of_true[[set_2_num]][[1]])),"_true_neg"=length(E(list_of_true[[set_2_num]][[2]]))))
          jacc <-common/(length(E(g_pos))+length(E(g_neg))+length(E(list_of_true[[set_2_num]][[1]]))+length(E(list_of_true[[set_2_num]][[2]]))-common)
          ran_res <- append(ran_res,round(jacc, digits=3))}
        res <- append(res,ran_res)}
    }}
  return (data.frame("bact"=bact, "jacc"=res))
}


Jaccard_rnd_network <- jaccard_rnd_vs_real_networks(list(WT_shuffle_networks,unevo_shuffle_networks, evo_shuffle_networks,evoWT_shuffle_networks,mut_shuffle_networks),list(WT,unevo,evo,evoWT,mut))

write_csv(Jaccard_rnd_network,"Jaccard_rnd_networks_pos_and_neg_PCLRC.csv")


