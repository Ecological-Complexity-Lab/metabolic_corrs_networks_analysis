#functions:

#organizing the data:

flattenCorrMatrix <- function(cormat, pmat,nmat) {
  ut <- upper.tri(cormat) 
  data.frame(
    M1 = rownames(cormat)[row(cormat)[ut]],
    M2 = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut],
    p_adj= p.adjust(pmat[ut], method = "fdr", n = length(pmat[ut])),
    n=nmat[ut]
  )}

#filter by n and p
filtering_corr <- function(results_cor, bact_name){
  df1 <- flattenCorrMatrix(results_cor$r, results_cor$P, results_cor$n)
  df_positive_p<-df1 %>% filter(cor > 0) 
  df_negative_p<-df1 %>% filter(cor < 0) 
  pos_and_neg_p_1 <- mapply(cbind, list(df1,df_positive_p,df_negative_p), "experiment"=bact_name, SIMPLIFY=F)
  
  return(pos_and_neg_p_1)
}

#get common_corrs with PCLRC corrs
common_corrs_PCLRC <- function(df1,df2){
  common_1= bind_rows(list(df1, df2)) %>%
    group_by(M1, M2) %>%filter(n()>1) %>%
    reframe(cor= cor)%>%drop_na()
  df_2_oppo <- setNames(df2,  c("M2", "M1"))
  common_2= bind_rows(list(df1, df_2_oppo)) %>%
    group_by(M1, M2) %>%filter(n()>1) %>%
    reframe(cor=cor)%>%drop_na()
  df_3=rbind(common_1,common_2)
  return(df_3)}

#common corrs between bacteria
common_corrs <- function(df1,df2){
  common_1= bind_rows(list(df1, df2)) %>%
    group_by(M1, M2) %>%filter(n()>1) %>%
    summarise(cor= sum(cor)/2)
  df_2_oppo <- setNames(df2,  c("M2", "M1", "cor"))
  common_2= bind_rows(list(df1, df_2_oppo)) %>%
    group_by(M1, M2) %>%filter(n()>1) %>%
    summarise(cor= sum(cor)/2)
  df_3=rbind(common_1,common_2)
  return(df_3)}


create_igraph <- function(df){
  before_igraph <- df%>% dplyr::select(from=M1, to=M2, weight=cor)
  igraph <-igraph::graph.data.frame(before_igraph ,directed=F) }


#jaccard between nets:
common_corrs_with_sign <- function(df1,df2){
  df1$sign=sign(df1$cor)
  df2$sign=sign(df2$cor)
  common_1= bind_rows(list(df1, df2)) %>%
    group_by(M1, M2, sign) %>%filter(n()>1) %>%
    summarise(cor= sum(cor)/2)
  df_2_oppo <- setNames(df2,  c("M2", "M1", "cor"))
  common_2= bind_rows(list(df1, df_2_oppo)) %>%
    group_by(M1, M2,sign) %>%filter(n()>1) %>%
    summarise(cor= sum(cor)/2)
  df_3=rbind(common_1,common_2)
  return(df_3)}


jaccard_bet_networks <- function(list_df){
  res <- as.data.frame(matrix(NA, nrow=length(list_df), ncol=length(list_df)))
  for (set_1_num in c(1:length(list_df))){
    for (set_2_num in c(1:length(list_df))){
      if (set_1_num!=set_2_num){
        common= common_corrs_with_sign(list_df[[set_1_num]],list_df[[set_2_num]])
        jacc <-nrow(common)/(nrow(list_df[[set_1_num]])+nrow(list_df[[set_2_num]])-nrow(common))
        res[set_1_num,set_2_num] <- round(jacc, digits=3)
      }}}
  return (res)
}

combined_column <- function(df){
  df2=df
  df2$combined <- paste(df2$M1, df2$M2, sep="_")
  return(df2$combined)
}


#add nodes that are singletons in one of the networks of data-set 1 or 2 the common network
add_singltons <- function(igraph_common,igraph1,igraph2){
  re <- intersect(names(V(igraph1)),names(V(igraph2)))
  print(length(re))
  rest <- setdiff(re,names(V(igraph_common)))
  print(length(rest))
  common_with_single<- add.vertices(igraph_common, length(rest), attr = list(name = rest))
  return(common_with_single)
}

#find common correlations between networks of different bacteria
core_between_networks <- function(list_df){
  n=1
  while(n+1<=length(list_df)){
    common= common_corrs(list_df[[n]],list_df[[n+1]])
    if (nrow(common)!=0){
      list_df[[n+1]]=common
      n=n+1
    }else{
      break}
  }
  print(n)
  return(common)
}

#shuffle the networks
set.seed(12)
shuffle_fun <- function(igraph){
  g1=rewire(igraph, with = keeping_degseq(niter = ecount(igraph)*10))
  E(g1)$weight <- sample(E(igraph)$weight,replace = FALSE)
  return(g1)
}

#compare jacc to random- used in a file called "pos_and_neg" in the folder "jacc_rnd_vs_true"
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

#get the results of the comparison of the true jaccard values to the random- #change the output to get the z values, p values or plots. 
z_score_of_jaccard_rnd <- function(jaccard_rnd_df, jaccard_true_df){
  res <- as.data.frame(matrix(NA, nrow=nrow(jaccard_true_df), ncol=ncol(jaccard_true_df)))
  plot_list=list()
  for (couple in c(1:nrow(jaccard_rnd_df))){
    print(couple)
    j_rnd_vector <- as.numeric(unlist(strsplit(jaccard_rnd_df[couple,2], ", ")))
    mean <- mean(j_rnd_vector)
    print(paste("mean_jacc=",mean))
    sd <- sd(j_rnd_vector)
    bact_couple <- unlist(strsplit(jaccard_rnd_df[couple,1], "_"))
    rep_str = c('1'='WT','2'='unevolved','3'='evolved', '4'='evolved-WT', 'rnd'='R','true'='T')
    bact_couple_names <- str_replace_all(jaccard_rnd_df[couple,1], rep_str)    
    bact_true <- as.numeric(bact_couple[4])
    print(bact_true)
    bact_rnd <- as.numeric(bact_couple[2])
    print(bact_rnd)
    jacc_true <- jaccard_true_df[bact_true,bact_rnd]
    print(jacc_true)
    z_val <- (jacc_true-mean)/sd
    p_val <-pnorm(z_val, mean = mean, sd = sd, lower.tail = FALSE)
    res[bact_true,bact_rnd]<- round(p_val, digits=3)
    d <- as.data.frame(j_rnd_vector)
    plot <- d%>% ggplot(aes(x=d[,1])) +geom_density( fill="dodgerblue", alpha=0.5)+geom_vline(xintercept=jacc_true, linewidth=1.5, color="red")+xlab(NULL)+ylab(NULL)+ggtitle(bact_couple_names)+theme_classic()+ theme(text = element_text(size = 15))
    
    plot_list <- append(plot_list,list(plot))
    
  }
  return(plot_list)
}


#nodes' roles'- by klil
library(magrittr)

z_score_multi <- function(i, modules, stats) {
  m_i <- modules$module[modules$state_node == i]
  k_m_i <- modules$k_m[modules$state_node == i]
  avg_k <- stats$k_m_avg[stats$module == m_i]
  sd_k <- stats$k_m_sd[stats$module == m_i]
  return ((k_m_i - avg_k)/sd_k)
}

c_score_multi <- function(i, modules, edges) {
  c_sum_i <- c()
  if (i %in% edges$sn_from) {
    i_edges <- filter(edges, sn_from == i)
    for (t in unique(i_edges$m_to)) {
      k_i_t <- sum(i_edges$m_to == t)
      k_i <- as.numeric(modules$k_total[modules$state_node == i])
      c_sum_i <- c(c_sum_i, (k_i_t / k_i)**2)
    }
  }
  if (i %in% edges$sn_to) {
    i_edges <- filter(edges, sn_to == i)
    for (t in unique(i_edges$m_from)) {
      k_i_t <- sum(i_edges$m_from == t)
      k_i <- as.numeric(modules$k_total[modules$state_node == i])
      c_sum_i <- c(c_sum_i, (k_i_t / k_i)**2)
    }
  }
  return(c_score_i <- as.numeric(1-sum(c_sum_i)))
}

assign_roles <-  function(modules) {
  output <- mutate(modules, 
                   role = case_when(modules$z_score <= 2.5 & modules$c_score <= 0.62 ~ 'peripheral',
                                    modules$z_score <= 2.5 & modules$c_score > 0.62 ~ 'connector',
                                    modules$z_score > 2.5 & modules$c_score <= 0.62 ~ 'module hub',
                                    modules$z_score > 2.5 & modules$c_score > 0.62 ~ 'network hub',
                                    is.nan(modules$z_score) == T & modules$c_score <= 0.62 ~ 'peripheral',
                                    is.nan(modules$z_score) == T & modules$c_score > 0.62 ~ 'connector',
                                    is.na(modules$z_score) == T & modules$c_score <= 0.62 ~ 'peripheral',
                                    is.na(modules$z_score) ==T & modules$c_score > 0.62 ~ 'connector'))
  return(output)
}


#first- create a df of the links with the module of each node- by state node
roles_fun <- function(modules_list,multi_net){
  expanded_intra <-multi_net$extended_ids %>% 
    left_join(modules_list, by = c('node_from' = 'node_id', 'layer_from' = 'layer_id')) %>% 
    left_join(modules_list, by = c('node_to' = 'node_id', 'layer_from' = 'layer_id')) %>% 
    dplyr::select(sn_from = state_node.x, sn_to = state_node.y, node_from, node_to, 
                  m_from = module.x, m_to = module.y, weight)
  main_graph <- graph.data.frame(expanded_intra[,c(1:2)], directed = F, vertices = modules_list$state_node) #v are state nodes
  
  modules_list %<>% mutate(k_total = igraph::degree(main_graph)) 
  expanded_list <- expanded_intra
  
  # degree within modules
  k_m <- data.frame(matrix(nrow = 0, ncol = 2))
  names(k_m) <- c('state_node', 'k_m')
  for (m in unique(modules_list$module)) {
    m_state_nodes <- filter(modules_list, module == m) 
    # edges where both nodes are in m 
    m_edges <- filter(expanded_list, m_from == m & m_to == m) 
    m_graph <- graph.data.frame(m_edges[,c(1:2)], directed = F, vertices = m_state_nodes$state_node)
    #V(m_graph)$type <- V(m_graph)$name %in% expanded_intra[,"sn_from"]
    # m_matrix <- as_incidence_matrix(m_graph, names = T, attr = 'weight', sparse = F)
    m_degrees <- data.frame(state_node = m_state_nodes$state_node, k_m = igraph::degree(m_graph)) 
    k_m <- rbind(k_m, m_degrees)
  }
  modules_list %<>% inner_join(k_m, by = "state_node") 
  
  
  
  # average degree in each module
  k_m_avg <- data.frame(matrix(nrow = 0, ncol = 3))
  for (m in unique(modules_list$module)) {
    k_m_avg %<>% rbind(c(m, mean(modules_list$k_m[modules_list$module == m]),
                         sd(modules_list$k_m[modules_list$module == m])))  
  }
  names(k_m_avg) <- c('module', 'k_m_avg', 'k_m_sd')
  
  
  modules_list %<>% arrange(state_node)
  Z <- lapply(modules_list$state_node, modules_list, k_m_avg, FUN = z_score_multi) 
  modules_list %<>% mutate(z_score = unlist(Z))
  
  C <- lapply(sort(modules_list$state_node), modules_list, expanded_intra, FUN = c_score_multi)
  modules_list %<>% mutate(c_score = unlist(C)) 
  
  super_list <- assign_roles(modules_list)
  write.csv(super_list, paste0("modularity/super_list_",str_split_1(deparse(substitute(multi_net)), "_")[1],"_PCLRC.csv"))
  
  #NOTE THAT SMALL MODULES THAT ARE DISCONNECTED FROM THE NETWORK CAN LEAD TO Z VALUE OF NA since the sd is 0- are not in the plot
  
  return(super_list)} 


create_NMI_matrix <- function(network_par_1,network_par_2){
  network_par_1$node_layer <- paste0(network_par_1$node_id,"_",network_par_1$layer_id)
  network_par_2$node_layer <- paste0(network_par_2$node_id,"_",network_par_2$layer_id)
  
  nodes_in_both <- intersect(network_par_1$node_layer,network_par_2$node_layer)
  network_par_1 <- network_par_1%>% drop_na()%>% filter( node_layer %in% nodes_in_both )
  print(nrow(network_par_1))
  network_par_2 <- network_par_2%>% drop_na()%>% filter( node_layer %in% nodes_in_both )
  print(nrow(network_par_2))
  network_par_1_modules <- network_par_1 %>% group_by(module) %>% count()
  network_par_2_modules <- network_par_2 %>% group_by(module) %>% count()
  mat <- matrix(0,nrow=max(as.numeric(network_par_1_modules$module)),ncol=max(as.numeric(network_par_2_modules$module)))
  for (metabolite in network_par_1$node_layer){ 
    mod_in_1 <- as.numeric(filter(network_par_1, node_layer==metabolite)$module)
    mod_in_2 <- as.numeric(filter(network_par_2, node_layer==metabolite)$module)
    mat[mod_in_1,mod_in_2]=mat[mod_in_1,mod_in_2]+1
  }
  return(mat)}

NMI_networks <- function(list_df){
  res <- as.data.frame(matrix(NA, nrow=length(list_df), ncol=length(list_df)))
  for (set_1_num in c(1:length(list_df))){
    for (set_2_num in c(1:length(list_df))){
      if (set_1_num!=set_2_num){#can change to > since the matrix is symmetric, but need all to compare to random
        NMI_mat <- create_NMI_matrix(list_df[[set_1_num]],list_df[[set_2_num]])
        delete_rows=c();delete_cols=c()
        NMI_mat=as.data.frame(NMI_mat) #now- delete rows/columns with sum==0 (since we deleted un-shared nodes and it couased the function NMI to stop)
        for (row in c(1:nrow(NMI_mat))){
          if (sum(NMI_mat[row,])==0){
            delete_rows <- append(delete_rows,row)
          }
        }
        for (col in c(1:ncol(NMI_mat))){
          if (sum(NMI_mat[,col])==0){
            delete_cols <- append(delete_cols,col)
          }
        }
        if (length(delete_rows)>0){
          NMI_mat <- NMI_mat[-delete_rows,]
        }
        if (length(delete_cols)>0){
          NMI_mat <- NMI_mat[,-delete_cols]
        }
        NMI_value <- infomapecology::NMI(NMI_mat)
        res[set_1_num,set_2_num] <- NMI_value
      }}}
  return (res)
}


node_info <- function(igraph){
  eigen=eigen_centrality(igraph, directed = FALSE, scale = FALSE,weights = abs(E(igraph)$weight), options = arpack_defaults)
  df2 <-data.frame("CC"= c(transitivity(igraph,type = "local",vids = V(igraph),weights = abs(E(igraph)$weight),isolates="zero")),
                   "bet"=c(betweenness(igraph,v = V(igraph),directed = FALSE, weights = 1/abs(E(igraph)$weight),normalized =F)),
                   "degree"=c(degree(igraph,v = V(igraph),mode ="all",loops = FALSE,normalized = F)),
                   "strength"=c(strength(igraph,vids = V(igraph),loops = FALSE)),
                   "eigen"=eigen$vector,
                   "nei"=influential::neighborhood.connectivity( igraph,vertices = V(igraph),mode = "all",verbose = FALSE))
  return(df2)}

MyMerge       <- function(x, y){
  df            <- merge(x, y, by= "row.names", all.x= T, all.y= T)
  rownames(df)  <- df$Row.names
  df$Row.names  <- NULL
  return(df)
}


CC_pos_all_high_degree <- CC_pos_all
for (met in c(1:nrow(CC_pos_all_high_degree))){
  for (bac in c(1:ncol(CC_pos_all_high_degree))){
    if (is.na(CC_pos_all_high_degree[met,bac])==FALSE & deg_pos_all[met,bac]<10){
      CC_pos_all_high_degree[met,bac] <- NA
    }}}
