# Load required libraries (used versions mentioned)

library(tidyverse) #2.0.0
library(readxl) #1.4.3
library(Hmisc) #4.7-1
library(ggplot2) #3.5.1
library(igraph) #1.5.1
library(data.table) #1.14.2
library(reshape2) #1.4.4
library(dplyr) #1.1.4
library(influential) #2.2.9
library(emln) #0.0.1
library(infomapecology)#2.0.1
library(magrittr) #2.0.3
library(combinat) #0.0-8
library(khroma) #1.12.0

# ----------------------------------------
# Step 1: Data Import and Pre-processing
# ----------------------------------------

# Define the experimental groups and batch names
groups <- c("WT", "unevo", "evolved", "evoWT", "mut")
batches <- c("585_1", "585_2")
sheet_names <- batches
data <- "Data/LC-MS_results.xlsx"

# Function to preprocess each sheet
process_data <- function(sheet_name, data_file = data) {
  df <- read_excel(data_file, sheet = sheet_name)  # Load Excel sheet into DataFrame
  
  # Build regex pattern for selecting group or 'cluster' (metabolic feature) columns
  pattern <- paste0("^(", paste(groups, collapse = "|"), "|cluster)")
  sample_cols <- grep(pattern, colnames(df), ignore.case = TRUE)
  
  # Replace zeros with NA in sample columns
  df[sample_cols][df[sample_cols] == 0] <- NA
  
  # Loop over each group and extract corresponding matrix
  lapply(groups, function(grp, df) {
    # Select all columns starting with the group name
    col_range <- df %>% select(starts_with(grp))
    
    # Add 'cluster' column before the group columns
    tbl <- col_range %>% add_column(cluster = df$cluster, .before = TRUE)
    
    # Set rownames to be the cluster identifiers
    rownames(tbl) <- tbl$cluster
    
    # Convert to matrix and apply log-transformation after transposition
    mat <- tbl %>%
      data.matrix() %>%
      t() %>%
      log()
    
    return(mat[-1,]) # Remove the 'cluster' row (which was transposed into a column)
  }, df) %>%
    setNames(groups)  # Assign group names to resulting list
}

# Apply the data processing to each batch (Excel sheet)
batch_data <- lapply(sheet_names, process_data)
names(batch_data) <- batches  # Name each processed batch by its original name

# ----------------------------------------
# Step 2: Correlation Analysis
# ----------------------------------------

# Step 2.1 - Define functions

# Function to compute Pearson's correlation values
compute_correlation <- function(data) {
  Hmisc::rcorr(as.matrix(data), type = "pearson")
}

# Function to convert upper triangle of correlation, p-value, and sample size matrices into a long-format dataframe
flatten_corr_matrix <- function(cormat, pmat, nmat) {
  upper_indices <- upper.tri(cormat)
  df <- data.frame(
    Feature1 = rownames(cormat)[row(cormat)[upper_indices]],
    Feature2 = rownames(cormat)[col(cormat)[upper_indices]],
    Correlation = cormat[upper_indices],
    P_value = pmat[upper_indices],
    Sample_Size = nmat[upper_indices]
  )
  # Adjust p-values for multiple testing using FDR
  df$P_value_adj <- p.adjust(df$P_value, method = "fdr", n = length(df$P_value))
  return(df)
}

# Function to filter correlations by sign (positive or negative)
filter_correlations <- function(cor_results, experiment_name) {
  df <- flatten_corr_matrix(cor_results$r, cor_results$P, cor_results$n)
  df <- df %>% mutate(Experiment = experiment_name)
  list(
    all = df,
    positive = df %>% filter(Correlation > 0),
    negative = df %>% filter(Correlation < 0)
  )
}

# Function to identify overlapping edges between standard correlation and PCLRC-filtered correlation edges
filter_pclrc <- function(cor_df, pclrc_df) {
  # Direct matches
  common_direct <- bind_rows(cor_df, pclrc_df) %>%
    group_by(Feature1, Feature2) %>%
    filter(n() > 1) %>%
    summarise(Correlation = first(Correlation), .groups = "drop") %>%
    drop_na()
  # Reverse matches
  pclrc_df_reversed <- pclrc_df %>% rename(Feature1 = Feature2, Feature2 = Feature1)
  common_reverse <- bind_rows(cor_df, pclrc_df_reversed) %>%
    group_by(Feature1, Feature2) %>%
    filter(n() > 1) %>%
    summarise(Correlation = first(Correlation), .groups = "drop") %>%
    drop_na()
  
  # Combine both direct and reverse matches
  bind_rows(common_direct, common_reverse)
}

# Step 2.2 - Correlation analysis for each batch using PCLRC-filtered edges

# NOTE: csv files of the PCLRC-filtered edges were generated with the R script from Di Cesare et al. (2022), GeroScience, 44(2), 1109???1128.
# This script retains only the edges that remain after applying the PCLRC filtering method.

# The function processes all datasets and compares them to precomputed PCLRC edge lists
perform_batch_correlation_analysis <- function(datasets, pclrc_edge_paths) {
  
  # Compute Pearson correlation matrices for each group dataset
  correlation_results <- lapply(datasets, compute_correlation)
  
  # Filter correlations by sign
  filtered_results <- mapply(filter_correlations, correlation_results, names(correlation_results), SIMPLIFY = FALSE)
  
  # Read PCLRC-filtered edges from CSVs for each group
  pclrc_edges <- lapply(pclrc_edge_paths, function(path) {
    readr::read_csv(path, show_col_types = FALSE) %>%
      select(2:3) %>%
      rename(Feature1 = 1, Feature2 = 2)
  })
  
  # Match positive correlations with PCLRC edges
  positive_overlap <- mapply(filter_pclrc,
                             lapply(filtered_results, `[[`, "positive"),
                             pclrc_edges,
                             SIMPLIFY = FALSE)
  
  # Match negative correlations with PCLRC edges
  negative_overlap <- mapply(filter_pclrc,
                             lapply(filtered_results, `[[`, "negative"),
                             pclrc_edges,
                             SIMPLIFY = FALSE)
  
  list(
    positive_overlap = positive_overlap,
    negative_overlap = negative_overlap,
    all_filtered_results = filtered_results
  )
}

# Batch 1 analysis

pclrc_files_batch1 <- list(
  WT = "Data/PCLRC_filtered_edges/WT_PCLRC_edges_1.csv",
  unevo = "Data/PCLRC_filtered_edges/unevo_PCLRC_edges_1.csv",
  evo = "Data/PCLRC_filtered_edges/evo_PCLRC_edges_1.csv",
  evoWT = "Data/PCLRC_filtered_edges/evoWT_PCLRC_edges_1.csv",
  mut = "Data/PCLRC_filtered_edges/mut_PCLRC_edges_1.csv"
)

batch1_results <- perform_batch_correlation_analysis(batch_data$`585_1`, pclrc_files_batch1)

# Batch 2 analysis

pclrc_files_batch2 <- list(
  WT = "Data/PCLRC_filtered_edges/WT_PCLRC_edges_2.csv",
  unevo = "Data/PCLRC_filtered_edges/unevo_PCLRC_edges_2.csv",
  evo = "Data/PCLRC_filtered_edges/evo_PCLRC_edges_2.csv",
  evoWT = "Data/PCLRC_filtered_edges/evoWT_PCLRC_edges_2.csv",
  mut = "Data/PCLRC_filtered_edges/mut_PCLRC_edges_2.csv"
)

batch2_results <- perform_batch_correlation_analysis(batch_data$`585_2`, pclrc_files_batch2)

# ----------------------------------------
# Step 3: Extract Common Correlations Between Batches
# ----------------------------------------

# Function to find shared correlations between two data frames,
# keeping only those with the same sign and averaging their values
common_corrs <- function(df1, df2) {
  df1 <- df1 %>% mutate(sign = sign(Correlation))
  df2 <- df2 %>% mutate(sign = sign(Correlation))
  
  # Find direct matches (same Feature1 and Feature2)
  common_1 <- bind_rows(df1, df2) %>%
    group_by(Feature1, Feature2,sign) %>%
    filter(n() > 1) %>%
    summarise(Correlation = mean(Correlation), .groups = "drop")
  
  # Find reverse matches (Feature1 <-> Feature2)
  df2_flipped <- setNames(df2, c("Feature2", "Feature1", "Correlation"))
  common_2 <- bind_rows(df1, df2_flipped) %>%
    group_by(Feature1, Feature2,sign) %>%
    filter(n() > 1) %>%
    summarise(Correlation = mean(Correlation), .groups = "drop")
  
  # Combine direct and reverse matches
  bind_rows(common_1, common_2)
}

# Wrapper function to apply "common_corrs" function to all groups and both positive/negative networks
extract_common_corrs_all_groups <- function(d1_pos, d1_neg, d2_pos, d2_neg) {
  groups <- names(d1_pos)
  map(groups, function(g) {
    list(
      positive = common_corrs(d1_pos[[g]], d2_pos[[g]]), # Common positive correlations
      negative = common_corrs(d1_neg[[g]], d2_neg[[g]])  # Common negative correlations
    )
  }) %>% set_names(groups) # Return named list by group
}

# Apply to both batches and store all common correlations per group
common_corrs_all <- extract_common_corrs_all_groups(
  batch1_results$positive_overlap, batch1_results$negative_overlap,
  batch2_results$positive_overlap, batch2_results$negative_overlap
)


# ----------------------------------------
# Step 4: Network Comparison and Similarity Analysis (Figure 1E)
# ----------------------------------------

# Function to compute Jaccard similarity index between all bacterial networks
jaccard_between_networks <- function(list_df) {
  n <- length(list_df)
  res <- matrix(NA, nrow = n, ncol = n)
  rownames(res) <- colnames(res) <- names(list_df)
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j) {
        # Get common edges between networks i and j
        common <- common_corrs(list_df[[i]], list_df[[j]])
        # Calculate Jaccard index: |A ??? B| / |A ??? B|
        jacc <- nrow(common) / (nrow(list_df[[i]]) + nrow(list_df[[j]]) - nrow(common))
        res[i, j] <- round(jacc, 3)
      } else {
        res[i, j] <- 1 # similarity of a network with itself
      }
    }
  }
  as.data.frame(res)
}

# Function to combine positive and negative edges for each group
combine_pos_neg <- function(common_list) {
  map(common_list, ~ bind_rows(.x$positive, .x$negative))
}

# Combine edges for each group
common_combined <- combine_pos_neg(common_corrs_all)

# Compute pairwise Jaccard similarity between all groups
Jaccard_networks <- jaccard_between_networks(common_combined)

#Calculate Jaccard similairty index between networks of the same bacterial group in the two batches

# Combine positive and negative edges per group for each batch
combine_pos_neg_per_batch <- function(pos_list, neg_list) {
  mapply(function(p, n) rbind(p, n), pos_list, neg_list, SIMPLIFY = FALSE)
}

batch1_combined <- combine_pos_neg_per_batch(batch1_results$positive_overlap, batch1_results$negative_overlap)
batch2_combined <- combine_pos_neg_per_batch(batch2_results$positive_overlap, batch2_results$negative_overlap)

# Match names to be sure
names(batch1_combined) <- names(batch2_combined) <- names(batch1_results$positive_overlap)

# Now compare the same group between batch1 and batch2
# Create a list of both versions for each group
combined_by_group <- mapply(function(b1, b2) list(batch1 = b1, batch2 = b2),
                            batch1_combined, batch2_combined,
                            SIMPLIFY = FALSE)

# Now apply Jaccard on each group separately
jaccard_results_betwenn_batches <- lapply(combined_by_group, jaccard_between_networks)


# ----------------------------------------
# Step 5: Network Construction and Summary Statistics (Table 2)
# ----------------------------------------

# Function to create an igraph network object from a filtered correlation data frame, using the common correlations as weights
create_igraph <- function(df) {
  edges <- df %>%
    dplyr::select(from = Feature1, to = Feature2, weight = Correlation)
  igraph::graph_from_data_frame(edges, directed = FALSE)
}

# Create igraph objects for all groups
igraphs_all <- lapply(common_corrs_all, function(group_list) {
  lapply(group_list, create_igraph)
})

# Create igraph objects for Batch 1, separately for positive and negative correlation edges
igraphs_batch_1_pos <- lapply(batch1_results$positive_overlap, create_igraph)
igraphs_batch_1_neg <- lapply(batch1_results$negative_overlap, create_igraph)

# Create igraph objects for Batch 2, separately for positive and negative correlation edges
igraphs_batch_2_pos <- lapply(batch2_results$positive_overlap, create_igraph)
igraphs_batch_2_neg <- lapply(batch2_results$negative_overlap, create_igraph)


# Step 5.1: Add singleton nodes to recover full node sets (all metabolic features)

# Function to add singleton vertices (nodes without edges) to the common network,
# to match the union of nodes present in the two batch-specific networks
add_singltons <- function(igraph_common,igraph1,igraph2){
  
  # Find intersection of node names present in both batch networks
  re <- intersect(names(V(igraph1)),names(V(igraph2)))
  
  # Determine nodes present in intersection but missing from the common network
  rest <- setdiff(re,names(V(igraph_common)))
  
  # Add missing nodes as vertices to the common network, preserving their names
  common_with_single<- add.vertices(igraph_common, length(rest), attr = list(name = rest))
  return(common_with_single)
}

# Wrapper function to apply "add_singltons" function across all groups for a (for positive or negative data)
add_singltons_all_nested <- function(igraphs_all, batch1_list, batch2_list, polarity) {
  lapply(names(igraphs_all), function(group_name) {
    ig_common <- igraphs_all[[group_name]][[polarity]]
    # Use 'unevo' group networks from both batches as reference for node sets, as it contains all nodes including singletons
    ig1 <- batch1_list[["unevo"]] 
    ig2 <- batch2_list[["unevo"]]
    # Add missing singleton nodes to the common network
    result <- add_singltons(ig_common, ig1, ig2)
    return(result)
  }) %>%
    setNames(names(igraphs_all))
}

# Add singleton nodes for all groups in positive networks
igraphs_all_pos_s <- add_singltons_all_nested(
  igraphs_all,
  igraphs_batch_1_pos,
  igraphs_batch_2_pos,
  polarity = "positive"
)

# Add singleton nodes for all groups in negative networks
igraphs_all_neg_s <- add_singltons_all_nested(
  igraphs_all,
  igraphs_batch_1_neg,
  igraphs_batch_2_neg,
  polarity = "negative"
)

# Step 5.2: Compute global network summary properties:
#number of nodes, number of edges, density, average degree, average correlation weight, network diameter,
#global clustering coefficient, average shortest path and average neighborhood connectivity.


# Initialize an empty data frame to hold summary statistics for positive and negative networks with singleton nodes
amount_of_v_and_e_with_singltons <- data.frame(
  Group = character(),
  V_pos = integer(),
  E_pos = integer(),
  dens_pos = numeric(),
  V_neg = integer(),
  E_neg = integer(),
  dens_neg = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each group and calculate number of vertices, edges, and edge density for positive and negative networks
for (g in groups) {
  pos_graph <- igraphs_all_pos_s[[g]]
  neg_graph <- igraphs_all_neg_s[[g]]
  
  amount_of_v_and_e_with_singltons <- rbind(
    amount_of_v_and_e_with_singltons,
    data.frame(
      Group = g,
      V_pos = length(V(pos_graph)),             # Number of nodes in positive network
      E_pos = length(E(pos_graph)),             # Number of edges in positive network
      dens_pos = round(edge_density(pos_graph), 4),  # Edge density of positive network
      V_neg = length(V(neg_graph)),             # Number of nodes in negative network
      E_neg = length(E(neg_graph)),             # Number of edges in negative network
      dens_neg = round(edge_density(neg_graph), 4)   # Edge density of negative network
    )
  )
}

# Set row names as group names and remove redundant Group column
rownames(amount_of_v_and_e_with_singltons) <- amount_of_v_and_e_with_singltons$Group
amount_of_v_and_e_with_singltons$Group <- NULL

# Repeat network summary statistics calculation WITHOUT singleton nodes included

# Initialize empty data frame for summary statistics without singletons
amount_of_v_and_e_without_singltons <- data.frame(
  Group = character(),
  V_pos = integer(),
  E_pos = integer(),
  dens_pos = numeric(),
  V_neg = integer(),
  E_neg = integer(),
  dens_neg = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each group and calculate stats for igraph objects without added singletons
for (g in groups) {
  pos_graph <- igraphs_all[[g]]$positive
  neg_graph <- igraphs_all[[g]]$negative
  
  amount_of_v_and_e_without_singltons <- rbind(
    amount_of_v_and_e_without_singltons,
    data.frame(
      Group = g,
      V_pos = length(V(pos_graph)),
      E_pos = length(E(pos_graph)),
      dens_pos = round(edge_density(pos_graph), 4),
      V_neg = length(V(neg_graph)),
      E_neg = length(E(neg_graph)),
      dens_neg = round(edge_density(neg_graph), 4)
    )
  )
}

# Set row names and remove Group column for the data frame without singletons
rownames(amount_of_v_and_e_without_singltons) <- amount_of_v_and_e_without_singltons$Group
amount_of_v_and_e_without_singltons$Group <- NULL


# Initialize empty data frames to store mean shortest path lengths
sp_avg_pos <- data.frame(bact = character(), mean_distance = numeric())
sp_avg_neg <- data.frame(bact = character(), mean_distance = numeric())

# Initialize empty data frames to store network diameters
diameter_pos <- data.frame(bact = character(), diameter = numeric())
diameter_neg <- data.frame(bact = character(), diameter = numeric())

# Initialize empty data frame to store global clustering coefficients
cc_df <- data.frame(bact = character(), CC_global_pos = numeric(), CC_global_neg = numeric())

# Initialize empty data frames to store average node degrees
avg_degree_pos <- data.frame(bact = character(), avg_degree = numeric())
avg_degree_neg <- data.frame(bact = character(), avg_degree = numeric())

# Initialize empty data frames to store average edge weights (correlation strengths)
avg_weight_pos <- data.frame(bact = character(), avg_weight = numeric())
avg_weight_neg <- data.frame(bact = character(), avg_weight = numeric())

# Initialize empty data frames to store average neighborhood connectivity
avg_nnc_pos <- data.frame(bact = character(), avg_nnc = numeric())
avg_nnc_neg <- data.frame(bact = character(), avg_nnc = numeric())


for (bact in names(igraphs_all)) {
  ig_pos <- igraphs_all[[bact]][[1]]  # Positive network igraph object
  ig_neg <- igraphs_all[[bact]][[2]]  # Negative network igraph object
  
  # Mean shortest path lengths
  sp_avg_pos <- rbind(sp_avg_pos, data.frame(
    bact = bact,
    mean_distance = mean_distance(ig_pos,
                                  weights = 1/E(ig_pos)$weight,
                                  directed = FALSE,
                                  unconnected = TRUE)
  ))
  
  sp_avg_neg <- rbind(sp_avg_neg, data.frame(
    bact = bact,
    mean_distance = mean_distance(ig_neg,
                                  weights = abs(1/E(ig_neg)$weight),
                                  directed = FALSE,
                                  unconnected = TRUE)
  ))
  
  # Network diameters
  diameter_pos <- rbind(diameter_pos, data.frame(
    bact = bact,
    diameter = diameter(ig_pos,
                        weights = 1/E(ig_pos)$weight,
                        directed = FALSE,
                        unconnected = TRUE)
  ))
  
  diameter_neg <- rbind(diameter_neg, data.frame(
    bact = bact,
    diameter = diameter(ig_neg,
                        weights = 1/abs(E(ig_neg)$weight),
                        directed = FALSE,
                        unconnected = TRUE)
  ))
  
  # Global clustering coefficient
  cc_pos <- transitivity(ig_pos, type = "global", vids = V(ig_pos),
                         weights = 1/abs(E(ig_pos)$weight), isolates = "zero")
  cc_neg <- transitivity(ig_neg, type = "global", vids = V(ig_neg),
                         weights = 1/abs(E(ig_neg)$weight), isolates = "zero")
  
  cc_df <- rbind(cc_df, data.frame(
    bact = bact,
    CC_global_pos = cc_pos,
    CC_global_neg = cc_neg
  ))
  
  # Average degree
  avg_degree_pos <- rbind(avg_degree_pos, data.frame(
    bact = bact,
    avg_degree = mean(degree(ig_pos))
  ))
  avg_degree_neg <- rbind(avg_degree_neg, data.frame(
    bact = bact,
    avg_degree = mean(degree(ig_neg))
  ))
  
  # Average edge weight (correlation strength)
  avg_weight_pos <- rbind(avg_weight_pos, data.frame(
    bact = bact,
    avg_weight = mean(E(ig_pos)$weight)
  ))
  avg_weight_neg <- rbind(avg_weight_neg, data.frame(
    bact = bact,
    avg_weight = mean(E(ig_neg)$weight)
  ))
  
  # Average neighborhood connectivity
  avg_nnc_pos <- rbind(avg_nnc_pos, data.frame(
    bact = bact,
    avg_nnc = mean(influential::neighborhood.connectivity(ig_pos))
  ))
  avg_nnc_neg <- rbind(avg_nnc_neg, data.frame(
    bact = bact,
    avg_nnc = mean(influential::neighborhood.connectivity(ig_neg))
  ))
  
}


# ----------------------------------------
# Step 6: Generate Shuffled Networks and Compare the Jaccard Indexes (Fig. 5S)
# ----------------------------------------


set.seed(12) # Set seed for reproducibility

# Function to shuffle edges of a network while preserving degree sequence,
# then randomly reassign edge weights by shuffling original weights
shuffle_fun <- function(igraph) {
  g1 = rewire(igraph, with = keeping_degseq(niter = ecount(igraph) * 10))
  E(g1)$weight <- sample(E(igraph)$weight, replace = FALSE)
  return(g1)
}

# Function to combine negative and positive networks of each bacterial strain
process_network <- function(positive_igraph, negative_igraph) {
  combined_igraph = positive_igraph %u% negative_igraph
  # handling missing edge weights by replacing NA with zero
  E(combined_igraph)$weight_1[is.na(E(combined_igraph)$weight_1)] <- 0
  E(combined_igraph)$weight_2[is.na(E(combined_igraph)$weight_2)] <- 0
  # taking weights from positive and negative networks
  E(combined_igraph)$weight <- E(combined_igraph)$weight_1 + E(combined_igraph)$weight_2
  return(combined_igraph)
}

# Function to combine positive and negative networks, then generate multiple shuffled versions
process_and_shuffle <- function(pos_igraph, neg_igraph, num_shuffles) {
  combined_igraph <- process_network(pos_igraph, neg_igraph)
  shuffled_networks <- replicate(num_shuffles, shuffle_fun(combined_igraph), simplify = FALSE)
  return(shuffled_networks)
}


# Apply the process_and_shuffle function to all groups.
# Here 10 shuffles are performed per group for runtime reasons (originally 1000)
all_shuffled_networks <- mapply(
  FUN = process_and_shuffle,
  pos_igraph = igraphs_all_pos_s,
  neg_igraph = igraphs_all_neg_s,
  num_shuffles = 10,
  SIMPLIFY = FALSE
)

# Function to compare Jaccard similarity between random (shuffled) and empirical networks
# Loops over pairs of networks, skipping self-comparisons.
# For each pair, computes Jaccard index based on overlapping positive and negative edges
jaccard_rnd_vs_real_networks <- function(list_of_rnd, list_of_true){
  res <- c()
  bact <- c()
  for (set_2_num in c(1:length(list_of_true))){ # iterate over "true" (empirical) networks
    for (set_1_num in c(1:length(list_of_rnd))){ # iterate over shuffled networks
      if (set_2_num!=set_1_num){
        ran_res <- c()
        for (ran_net in c(1:length(list_of_rnd[[set_1_num]]))){ # iterate over shuffled replicates
          rnd=list_of_rnd[[set_1_num]][[ran_net]]
          
          # Separate positive and negative edges by weight sign
          g_pos <- delete.edges(rnd, which(E(rnd)$weight <0))
          g_neg <- delete.edges(rnd, which(E(rnd)$weight >0))
          
          # Count common edges between shuffled and true networks for positive and negative graphs
          pos_common=length(E(g_pos %s% list_of_true[[set_2_num]][[1]]))
          neg_common=length(E(g_neg %s% list_of_true[[set_2_num]][[2]]))
          common= pos_common+neg_common
          bact <- append(bact,paste0("rnd_",set_1_num,"_","true_",set_2_num))

          # Calculate Jaccard index for combined edges
          jacc <-common/(length(E(g_pos))+length(E(g_neg))+length(E(list_of_true[[set_2_num]][[1]]))+length(E(list_of_true[[set_2_num]][[2]]))-common)
          ran_res <- append(ran_res,round(jacc, digits=3))}
        res <- append(res,ran_res)}
    }}
  return (data.frame("bact"=bact, "jacc"=res))
}

# Combine positive and negative igraphs with singletons into list pairs for input to comparison
igraphs_all_s <- Map(function(pos, neg) list(pos, neg), igraphs_all_pos_s, igraphs_all_neg_s)

# Run Jaccard similarity comparison between shuffled and empirical networks
Jaccard_rnd_network <- jaccard_rnd_vs_real_networks(all_shuffled_networks, igraphs_all_s)

# Aggregate results by each pair of networks (random vs true)
jacc_rnd_pos_and_neg1 <- aggregate(jacc ~ bact, Jaccard_rnd_network, paste, collapse = ", ")


# Function to compute z-scores and p-values comparing empirical Jaccard indices
# against distributions of Jaccard indices from shuffled networks

z_score_of_jaccard_rnd <- function(jaccard_rnd_df, jaccard_true_df) {
  res_z <- matrix(NA, nrow = nrow(jaccard_true_df), ncol = ncol(jaccard_true_df),
                  dimnames = list(rownames(jaccard_true_df), colnames(jaccard_true_df)))
  res_p <- res_z
  plot_list <- list()
  
  # Mapping of string replacements for nicer plot titles
  rep_str <- c('1'='WT','2'='unevolved','3'='evolved',
               '4'='evolved-WT', '5'='mut', 'rnd'='R', 'true'='T')
  
  for (i in seq_len(nrow(jaccard_rnd_df))) {
    row_name <- jaccard_rnd_df[i, 1]
    j_rnd_vector <- as.numeric(unlist(strsplit(jaccard_rnd_df[i, 2], ", ")))
    
    m <- mean(j_rnd_vector)
    s <- sd(j_rnd_vector)
    
    # Extract indices from the string (format "rnd_X_true_Y")
    bact_couple <- unlist(strsplit(row_name, "_"))
    bact_true <- as.integer(bact_couple[4])
    bact_rnd <- as.integer(bact_couple[2])
    
    jacc_true <- jaccard_true_df[bact_true, bact_rnd]
    
    z_val <- (jacc_true - m) / s
    p_val <- pnorm(z_val, mean = m, sd = s, lower.tail = FALSE)
    
    res_z[bact_true, bact_rnd] <- round(z_val, 3)
    res_p[bact_true, bact_rnd] <- p_val
    
    # Create nicer label for plots by replacing codes
    bact_couple_names <- stringr::str_replace_all(row_name, rep_str)
    d <- data.frame(Jaccard = j_rnd_vector)
    
    plot <- ggplot(d, aes(x = Jaccard)) +
      geom_histogram(binwidth = 0.0005, fill = "dodgerblue", alpha = 0.5, boundary = 0) +
      geom_vline(xintercept = jacc_true, linewidth = 1.5, color = "red") +
      xlab(NULL) + ylab(NULL) + ggtitle(bact_couple_names) +
      theme_classic() + theme(text = element_text(size = 12))
    
    plot_list[[length(plot_list) + 1]] <- plot
  }
  
  return(list(plots = plot_list, z_scores = res_z, p_values = res_p))
}

# Run z-score and p-value calculations and generate density plots comparing shuffled vs empirical Jaccard indices
jacc_rnd_vs_true_results <- z_score_of_jaccard_rnd(jacc_rnd_pos_and_neg1,Jaccard_networks)

#Fig 5S
s <- arrangeGrob(grobs = jacc_rnd_vs_true_results[[1]])
figure <- annotate_figure(s, left = textGrob("Count", rot = 90, vjust = 1, gp = gpar(cex = 1.7)),
                          bottom = textGrob("Jaccard Values", gp = gpar(cex = 1.7)),fig.lab.size=0.5,fig.lab.face = "bold")

# ----------------------------------------
# Step 7: Create Multilayer Networks Using the "emln" Package
# ----------------------------------------

# Note: The emln package requires positive weights for proper functionality,
# so the absolute value of the negative edge weights is used.

# Initialize list to store multilayer networks for each group
multilayer_networks <- list()

# Loop over each bacterial group to build multilayer networks
for (group in names(igraphs_all_s)) {
  pos_graph <- igraphs_all_s[[group]][[1]]
  neg_graph <- igraphs_all_s[[group]][[2]]
  
  # Convert negative edge weights to positive by taking absolute value
  E(neg_graph)$weight <- abs(E(neg_graph)$weight)  # Ensure positive weights for negative layer
  
  # Extract adjacency matrices weighted by edge weights for both layers
  pos_mat <- as.matrix(pos_graph, attr = "weight", sparse = FALSE)
  neg_mat <- as.matrix(neg_graph, attr = "weight", sparse = FALSE)
  
  # Create multilayer network object from positive and negative adjacency matrices
  multilayer_networks[[group]] <- emln::create_multilayer_network(
    list_of_layers = list(pos_mat, neg_mat),bipartite = FALSE,directed = FALSE)
}

# ----------------------------------------
# Step 8: Modularity Using 'Infomap' Package
# ----------------------------------------

# This step performs community detection on multilayer networks using the Infomap algorithm.
# Note: This analysis is computationally intensive and is recommended to run on an HPC cluster.
# Ensure 'Infomap' package is installed and functional (use check_infomap() to verify).
# Optionally, you can skip running the analysis here and load previously saved results at the end of this section.

# Function to run Infomap community detection on a multilayer network object
modularity_using_infomap_multi <- function(multilayer_net) {
  # Create the multilayer object compatible with Infomap
  multi_obj <- infomapecology::create_multilayer_object(
    extended = multilayer_net$extended_ids,
    nodes = multilayer_net$nodes,
    layers = multilayer_net$layers,
    intra_output_extended = FALSE
  )
  multi_obj$inter <- NULL   # Remove inter-layer links if none exist (optional)
  
  #Run the Infomap multilayer algorithm with specified parameters:
    
    # - relax = TRUE: allows node state relaxations across layers
    # - flow_model = 'undirected': undirected flow assumed
    # - silent = TRUE: suppresses console output during run
    # - trials = 100: number of algorithm trials to optimize results
    # - seed = 1234: fixed seed for reproducibility
    # - temporal_network = FALSE: treat as static multilayer network
    
  multi_info <- run_infomap_multilayer(
    M = multi_obj,
    relax = TRUE,
    flow_model = 'undirected',
    silent = TRUE,
    trials = 100,
    seed = 1234,
    temporal_network = FALSE
  )
  
  # Process module assignments returned by Infomap
  module_list <- multi_info$modules %>%
    mutate(state_node = seq_len(n())) %>% # Add sequential node ID
    arrange(module) %>% # Sort by module number
    select(-starts_with("."))  # Remove columns starting with '.'
  
  # Return list containing:
  # - full Infomap output object
  # - multilayer object used
  # - cleaned modules data frame
  return(list(info = multi_info, obj = multi_obj, modules = module_list))
}

# Apply the modularity detection function to all multilayer networks
modularity_results <- lapply(multilayer_networks, modularity_using_infomap_multi)

# Extract separate components from the results for each bacterial group:

# 1) Infomap output objects ('info'):

WT_info_object_multi <- modularity_results[["WT"]]$info
unevo_info_object_multi <- modularity_results[["unevo"]]$info
evo_info_object_multi <- modularity_results[["evo"]]$info
evoWT_info_object_multi <- modularity_results[["evoWT"]]$info
mut_info_object_multi <- modularity_results[["mut"]]$info

# 2) Multilayer network objects ('obj') that were analyzed by Infomap

WT_multi_obj_multi <- modularity_results[["WT"]]$obj
unevo_multi_obj_multi <- modularity_results[["unevo"]]$obj
evo_multi_obj_multi <- modularity_results[["evo"]]$obj
evoWT_multi_obj_multi <- modularity_results[["evoWT"]]$obj
mut_multi_obj_multi <- modularity_results[["mut"]]$obj

# 3) Module assignment data frames for each group (note: singletons may have NA values)

all_modules <- lapply(modularity_results, function(x) x$modules)

WT_modules <- modularity_results[["WT"]]$modules
unevo_modules <- modularity_results[["unevo"]]$modules
evo <- modularity_results[["evo"]]$modules
evoWT_modules <- modularity_results[["evoWT"]]$modules
mut_modules <- modularity_results[["mut"]]$modules

# ----------------------------------------
# Optional: Load precomputed Infomap modularity results instead of running above (to save time)

# 1) Load Infomap output objects ('info')

WT_info_object_multi <- read_rds("Data/multi_networks_results/WT_info_object.RData")
unevo_info_object_multi <- read_rds("Data/multi_networks_results/unevo_info_object.RData")
evo_info_object_multi <- read_rds("Data/multi_networks_results/evo_info_object.RData")
evoWT_info_object_multi <- read_rds("Data/multi_networks_results/evoWT_info_object.RData")
mut_info_object_multi <- read_rds("Data/multi_networks_results/mut_info_object.RData")

# 2) Load multilayer network objects ('obj')

WT_multi_obj_multi <- read_rds("Data/multi_networks_results/WT_multi_object.RData")
unevo_multi_obj_multi <- read_rds("Data/multi_networks_results/evo_multi_object.RData")
evo_multi_obj_multi <- read_rds("Data/multi_networks_results/evo_multi_object.RData")
evoWT_multi_obj_multi <- read_rds("Data/multi_networks_results/evoWT_multi_object.RData")
mut_multi_obj_multi <- read_rds("Data/multi_networks_results/mut_multi_object.RData")

# 3) Load module assignment data frames (note removal of hidden columns starting with '.')

WT_modules <- read_csv("Data/multi_networks_results/WT_module_list.csv") %>%select(-starts_with("."))
unevo_modules <- read_csv("Data/multi_networks_results/unevo_module_list.csv")%>%select(-starts_with("."))
evo_modules <- read_csv("Data/multi_networks_results/evo_module_list.csv")%>%select(-starts_with("."))
evoWT_modules <- read_csv("Data/multi_networks_results/evoWT_module_list.csv")%>%select(-starts_with("."))
mut_modules <- read_csv("Data/multi_networks_results/mut_module_list.csv")%>%select(-starts_with("."))

# ----------------------------------------

# ----------------------------------------
# Step 9: Assigning Node Roles in Multilayer Networks
# ----------------------------------------

# Node roles classification is based on two key metrics from Guimer?? & Amaral (2005):
# 1. Within-module degree (Z): measures how well a node is connected to others in its own module.
# 2. Participation coefficient (C): measures how a node's connections are distributed among different modules.
# Nodes are categorized into roles such as peripheral, connector, module hub, or network hub based on thresholds of Z and C.

# Function to calculate the Z-score for a given node 'i'
# Z-score quantifies how the node's within-module degree compares to the average within-module degree in that module
z_score_multi <- function(i, modules, stats) {
  m_i <- modules$module[modules$state_node == i]
  k_m_i <- modules$k_m[modules$state_node == i]
  avg_k <- stats$k_m_avg[stats$module == m_i]
  sd_k <- stats$k_m_sd[stats$module == m_i]
  return ((k_m_i - avg_k)/sd_k)
}

# Function to calculate the C-score (participation coefficient) for a node 'i'
# Measures distribution of node's links across different modules
c_score_multi <- function(i, modules, edges) {
  c_sum_i <- c()
  
  # Calculate contribution from edges where node i is the source
  if (i %in% edges$sn_from) {
    i_edges <- filter(edges, sn_from == i)
    for (t in unique(i_edges$m_to)) {
      k_i_t <- sum(i_edges$m_to == t) # Number of edges from node i to module t
      k_i <- as.numeric(modules$k_total[modules$state_node == i]) # Total degree of node i
      c_sum_i <- c(c_sum_i, (k_i_t / k_i)**2)  # Sum squared proportion of links to module t
    }
  }
  # Calculate contribution from edges where node i is the target
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

# Function to assign categorical roles to nodes based on Z and C thresholds
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


# Function to build a data frame of edges with module assignments for the connected nodes
# Then calculates within-module degrees and assigns roles to nodes
roles_fun <- function(modules_list,multi_net){
  
  # Join module assignments to edges (state nodes and modules for source and target)
  expanded_intra <-multi_net$extended_ids %>% 
    left_join(modules_list, by = c('node_from' = 'node_id', 'layer_from' = 'layer_id')) %>% 
    left_join(modules_list, by = c('node_to' = 'node_id', 'layer_from' = 'layer_id')) %>% 
    dplyr::select(sn_from = state_node.x, sn_to = state_node.y, node_from, node_to, 
                  m_from = module.x, m_to = module.y, weight)
  
  # Create a graph of all state nodes to compute degrees
  main_graph <- graph.data.frame(expanded_intra[,c(1:2)], directed = F, vertices = modules_list$state_node) #v are state nodes
  
  # Calculate total degree per node
  modules_list %<>% mutate(k_total = igraph::degree(main_graph)) 
  expanded_list <- expanded_intra
  
  # Calculate within-module degree (k_m) for each node
  k_m <- data.frame(matrix(nrow = 0, ncol = 2))
  names(k_m) <- c('state_node', 'k_m')
  for (m in unique(modules_list$module)) {
    m_state_nodes <- filter(modules_list, module == m) 
    
    # Select edges where both nodes belong to the same module m
    m_edges <- filter(expanded_list, m_from == m & m_to == m) 
    
    # Create graph for module m
    m_graph <- graph.data.frame(m_edges[,c(1:2)], directed = F, vertices = m_state_nodes$state_node)
    
    # Calculate degrees within module m
    m_degrees <- data.frame(state_node = m_state_nodes$state_node, k_m = igraph::degree(m_graph)) 
    k_m <- rbind(k_m, m_degrees)
  }
  
  # Join within-module degree back to modules_list
  modules_list %<>% inner_join(k_m, by = "state_node") 
  
  
  # Compute average and SD of within-module degrees per module for Z-score calculation
  k_m_avg <- data.frame(matrix(nrow = 0, ncol = 3))
  for (m in unique(modules_list$module)) {
    k_m_avg %<>% rbind(c(m, mean(modules_list$k_m[modules_list$module == m]),
                         sd(modules_list$k_m[modules_list$module == m])))  
  }
  names(k_m_avg) <- c('module', 'k_m_avg', 'k_m_sd')
  
  
  modules_list %<>% arrange(state_node)
  
  # Calculate Z-scores for all nodes
  Z <- lapply(modules_list$state_node, modules_list, k_m_avg, FUN = z_score_multi) 
  modules_list %<>% mutate(z_score = unlist(Z))
  
  # Calculate C-scores for all nodes
  C <- lapply(sort(modules_list$state_node), modules_list, expanded_intra, FUN = c_score_multi)
  modules_list %<>% mutate(c_score = unlist(C)) 
  
  # Assign categorical roles based on Z and C scores
  super_list <- assign_roles(modules_list)
  
  return(super_list)} 

# Assign roles to nodes for each group multilayer network using the modules identified by Infomap

WT_roles <- roles_fun(WT_modules,multilayer_networks[["WT"]])
unevo_roles <- roles_fun(unevo_modules,multilayer_networks[["unevo"]])
evo_roles <- roles_fun(evo_modules,multilayer_networks[["evolved"]])
evoWT_roles <- roles_fun(evoWT_modules,multilayer_networks[["evoWT"]])
mut_roles <- roles_fun(mut_modules,multilayer_networks[["mut"]])


# ----------------------------------------
# Step 10: Normalized Mutual Information (NMI) Calculation (Fig. 3B)
# ----------------------------------------

# NMI quantifies similarity between two network partitions (module assignments).

# Function to create the confusion matrix between two network partitions (where entries count the number of shared nodes between modules).

create_NMI_matrix <- function(network_par_1,network_par_2){
  
  # Create unique identifiers combining node_id and layer_id for matching nodes
  network_par_1$node_layer <- paste0(network_par_1$node_id,"_",network_par_1$layer_id)
  network_par_2$node_layer <- paste0(network_par_2$node_id,"_",network_par_2$layer_id)
  
  # Identify nodes present in both partitions
  nodes_in_both <- intersect(network_par_1$node_layer,network_par_2$node_layer)
  
  # Filter to only include nodes present in both partitions and remove NA modules
  network_par_1 <- network_par_1%>% drop_na()%>% filter( node_layer %in% nodes_in_both )
  network_par_2 <- network_par_2%>% drop_na()%>% filter( node_layer %in% nodes_in_both )
  
  # Count modules in each partition to determine matrix size
  network_par_1_modules <- network_par_1 %>% group_by(module) %>% count()
  network_par_2_modules <- network_par_2 %>% group_by(module) %>% count()
  
  # Initialize confusion matrix with zeros
  mat <- matrix(0,nrow=max(as.numeric(network_par_1_modules$module)),ncol=max(as.numeric(network_par_2_modules$module)))
  
  # For each node, increment the cell corresponding to its module assignments in both partitions
  for (metabolite in network_par_1$node_layer){ 
    mod_in_1 <- as.numeric(filter(network_par_1, node_layer==metabolite)$module)
    mod_in_2 <- as.numeric(filter(network_par_2, node_layer==metabolite)$module)
    mat[mod_in_1,mod_in_2]=mat[mod_in_1,mod_in_2]+1
  }
  return(mat)}


# Function to compute pairwise NMI values for a list of network partitions

NMI_networks <- function(list_df){
  res <- as.data.frame(matrix(NA, nrow=length(list_df), ncol=length(list_df)))
  for (set_1_num in c(1:length(list_df))){
    for (set_2_num in c(1:length(list_df))){
      if (set_1_num!=set_2_num){#can change to > since the matrix is symmetric, but need all to compare to random
        
        # Compute confusion matrix between two partitions
        NMI_mat <- create_NMI_matrix(list_df[[set_1_num]],list_df[[set_2_num]])
        
        # Remove empty rows and columns (modules with zero shared nodes)
        delete_rows <- c()
        delete_cols <- c()
        
        NMI_mat=as.data.frame(NMI_mat)
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
        
        # Calculate NMI value using infomapecology package function
        NMI_value <- infomapecology::NMI(NMI_mat)
        res[set_1_num,set_2_num] <- NMI_value
      }}}
  return (res)
}

# Compute the NMI matrix comparing module assignments across all bacterial groups/networks
NMI_all_bac <- NMI_networks(list(WT_modules,unevo_modules,evo_modules,evoWT_modules, mut_modules))%>%`colnames<-`(c("WT","unevo","evo","evoWT","mut"))%>%`row.names<-`(c("WT","unevo","evo","evoWT","mut"))


# ----------------------------------------
# Step 11: Node Level Analysis - Node Parameters
# ----------------------------------------

# Function to calculate node-level metrics for a given igraph object:
# - Local clustering coefficient (CC)
# - Betweenness centrality (bet)
# - Degree (number of connections)
# - Neighborhood connectivity (average degree of neighbors)


node_info <- function(igraph){
  df2 <-data.frame("CC"= c(transitivity(igraph,type = "local",vids = V(igraph),weights = 1/abs(E(igraph)$weight),isolates="zero")),
                   "bet"=c(betweenness(igraph,v = V(igraph),directed = FALSE, weights = 1/abs(E(igraph)$weight),normalized =T)),
                   "degree"=c(igraph::degree(igraph,v = V(igraph),mode ="all",loops = FALSE,normalized = F)),
                   "nei"=influential::neighborhood.connectivity( igraph,vertices = V(igraph),mode = "all",verbose = F))
  return(df2)}

# Apply node_info to each group's positive and negative networks, nested in igraphs_all list
node_info_all <- lapply(igraphs_all, function(group) {
  lapply(group, node_info)
})


# Extract positive node metrics for each group into separate named lists

pos_node_infos <- list(
  WT = node_info_all$WT$positive,
  unevo = node_info_all$unevo$positive,
  evo = node_info_all$evolved$positive,
  evoWT = node_info_all$evoWT$positive,
  mut = node_info_all$mut$positive
)

# Extract negative node metrics for each group similarly

neg_node_infos <- list(
  WT = node_info_all$WT$negative,
  unevo = node_info_all$unevo$negative,
  evo = node_info_all$evolved$negative,
  evoWT = node_info_all$evoWT$negative,
  mut = node_info_all$mut$negative
)

# Helper function to merge multiple data frames by row names
# Ensures all row names are preserved and returned as rownames (not a column)

merge_by_rownames <- function(dfs) {
  Reduce(function(x, y) merge(x, y, by = "row.names", all = TRUE) %>% 
           { rownames(.) <- .$Row.names; .[ , -1, drop = FALSE] },
         dfs)
}

# Extract a specific metric column from each group's node info, then merge by rownames

extract_and_merge_metric <- function(igraphs_list, metric) {
  dfs <- lapply(igraphs_list, function(df) df[, metric, drop = FALSE])
  names(dfs) <- names(igraphs_list)
  result <- merge_by_rownames(dfs)
  colnames(result) <- names(igraphs_list)
  return(result)
}

# Extract and merge positive network metrics for all groups

deg_pos_all   <- extract_and_merge_metric(pos_node_infos, "degree")
CC_pos_all    <- extract_and_merge_metric(pos_node_infos, "CC")
bet_pos_all   <- extract_and_merge_metric(pos_node_infos, "bet")
neib_pos_all  <- extract_and_merge_metric(pos_node_infos, "nei")

# Extract and merge negative network metrics for all groups
deg_neg_all   <- extract_and_merge_metric(neg_node_infos, "degree")
CC_neg_all    <- extract_and_merge_metric(neg_node_infos, "CC")
bet_neg_all   <- extract_and_merge_metric(neg_node_infos, "bet")
neib_neg_all  <- extract_and_merge_metric(neg_node_infos, "nei")

