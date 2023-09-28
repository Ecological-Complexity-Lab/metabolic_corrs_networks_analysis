# 7. Species roles in the connected network:
# counting interlayer edges as part of degree
library(infomapecology)
library(igraph)
library(tidyverse)
library(magrittr)
#library(betalink)
library(readxl)
library(ggalluvial)
library(combinat)
library(ggpubr)
library(grid)

check_infomap()

#in the main script I created multi_nets using the function from the package 'emln'- 'create_multilayer_network'. Used those networks to create_multi_obj to run the run_infomap_multi on.
#using infomap- modularity

modularity_using_infomap_multi <- function(multilayer_net){
  multi_obj <- create_multilayer_object(extended = multilayer_net$extended_ids, #creating multilayer from package infomap_ecology to use in "run_infomap_multilayer"
                                        nodes = multilayer_net$nodes,
                                        layers=multilayer_net$layers,intra_output_extended = F)
  multi_obj$inter <- NULL
  
  multi_info <- run_infomap_multilayer(M = multi_obj,relax = T, 
                                       flow_model = 'undirected', silent = T, 
                                       trials = 100, seed = 1234, temporal_network = F)
  
  module_list <- multi_info$modules 
  module_list  <- module_list %>% 
    mutate(state_node = c(1:length(module_list$node_id))) %>% 
    arrange(module)
  
  write.csv(module_list, paste0('modularity/multi_nets/module_list_',deparse(substitute(multilayer_net)),'_PCLRC.csv'))
  return(list(multi_info,multi_obj))}


WT_multilayer_net <- read_rds("modularity/multi_nets/WT_multilayer_PCLRC.RData") #this is the multilayer from package emln
unevo_multilayer_net <- read_rds("modularity/multi_nets/unevo_multilayer_PCLRC.RData") #this is the multilayer from package emln
evo_multilayer_net <- read_rds("modularity/multi_nets/evo_multilayer_PCLRC.RData") #this is the multilayer from package emln
evoWT_multilayer_net <- read_rds("modularity/multi_nets/evoWT_multilayer_PCLRC.RData") #this is the multilayer from package emln
mut_multilayer_net <- read_rds("modularity/multi_nets/mut_multilayer_PCLRC.RData") #this is the multilayer from package emln

WT <- modularity_using_infomap_multi(WT_multilayer_net)
unevo <- modularity_using_infomap_multi(unevo_multilayer_net)
evo <- modularity_using_infomap_multi(evo_multilayer_net)
evoWT <- modularity_using_infomap_multi(evoWT_multilayer_net)
mut <- modularity_using_infomap_multi(mut_multilayer_net)


WT_info_object_multi <- WT[[1]]
unevo_info_object_multi <- unevo[[1]]
evo_info_object_multi <- evo[[1]]
evoWT_info_object_multi <- evoWT[[1]]
mut_info_object_multi <- mut[[1]]

saveRDS(WT_info_object_multi,"modularity/WT_info_object_multi_PCLRC.RData")
saveRDS(unevo_info_object_multi,"modularity/unevo_info_object_multi_PCLRC.RData")
saveRDS(evo_info_object_multi,"modularity/evo_info_object_multi_PCLRC.RData")
saveRDS(evoWT_info_object_multi,"modularity/evoWT_info_object_multi_PCLRC.RData")
saveRDS(mut_info_object_multi,"modularity/mut_info_object_multi_PCLRC.RData")

WT_multi_obj_multi <- WT[[2]]
unevo_multi_obj_multi <- unevo[[2]]
evo_multi_obj_multi <- evo[[2]]
evoWT_multi_obj_multi <- evoWT[[2]]
mut_multi_obj_multi <- mut[[2]]


saveRDS(WT_multi_obj_multi,"modularity/WT_multi_obj_multi_PCLRC.RData")
saveRDS(unevo_multi_obj_multi,"modularity/unevo_multi_obj_multi_PCLRC.RData")
saveRDS(evo_multi_obj_multi,"modularity/evo_multi_obj_multi_PCLRC.RData")
saveRDS(evoWT_multi_obj_multi,"modularity/evoWT_multi_obj_multi_PCLRC.RData")
saveRDS(mut_multi_obj_multi,"modularity/mut_multi_obj_multi_PCLRC.RData")


