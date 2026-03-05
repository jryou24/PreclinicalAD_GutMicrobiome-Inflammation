calc_secom_network_properties <- function(secom_result, overlap_threshold = 10, cor_threshold = 0.5, seed_value = 123){
  library(igraph)  
  set.seed(seed_value)
  corr_linear = secom_result$corr_fl
  cooccur_linear = secom_result$mat_cooccur
  
  # Filter by co-occurrence
  corr_linear[cooccur_linear < overlap_threshold] = 0
  
  df_linear = data.frame(get_upper_tri(corr_linear)) %>%
    rownames_to_column("var1") %>%
    pivot_longer(cols = -var1, names_to = "var2", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(value = round(value, 2))
  
  tax_name = sort(union(df_linear$var1, df_linear$var2))
  df_linear$var1 = factor(df_linear$var1, levels = tax_name)
  df_linear$var2 = factor(df_linear$var2, levels = tax_name)
  
  # Create secom adjacency matrix to igraph object
  adj_matrix <- ifelse(abs(corr_linear) > 0.5, abs(corr_linear), 0)
  
  # Convert adjacency matrix to an igraph object
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE, weighted = TRUE)
  
  # Remove isolated nodes (nodes with degree 0)
  g <- delete_vertices(g, V(g)[degree(g) == 0])
  n_nodes <- vcount(g)
  
  # Extract the edge list for gephi
  edge_list <- as.data.frame(as_edgelist(g))
  colnames(edge_list) <- c("Source", "Target")
  
  # Extract indices of edges from adjacency matrix
  edge_indices <- which(upper.tri(adj_matrix) & adj_matrix, arr.ind = TRUE)
  
  # Extract corresponding correlation values
  edge_weights <- corr_linear[edge_indices]
  
  # Create weight column in gephi edge_list
  edge_direction <- ifelse(edge_weights > 0, "positive", "negative")
  edge_list$relation <- edge_direction
  edge_list$Weight <- abs(edge_weights)
  
  # Assign absolute weights to the edges in the igraph object
  E(g)$weight <- abs(edge_weights)
  
  # Now Calculate centrality measures
  degree_centrality <- degree(g)
  deg_norm <- degree_centrality / (n_nodes - 1) # normalized degree 
  betweenness_centrality <- betweenness(g, weights = E(g)$weight)
  bet_norm <- betweenness_centrality / (((n_nodes - 1) * (n_nodes - 2)) / 2) # normalize by max possible
  closeness_centrality <- closeness(g, weights = E(g)$weight, normalized = TRUE)

  
  local_eff_manual <- function(g) {
    sapply(V(g), function(v) {
      nbs <- neighbors(g, v)
      if (length(nbs) < 2) return(0)  # no subgraph possible, efficiency = 0
      
      # Subgraph induced by neighbors
      subg <- induced_subgraph(g, nbs)
      
      # Compute distances among neighbors
      d <- distances(subg)
      inv_d <- 1 / d
      diag(inv_d) <- 0
      
      # Count valid paths (avoid Inf)
      valid_paths <- inv_d[is.finite(inv_d)]
      
      if (length(valid_paths) == 0) return(0)  # no connections among neighbors
      
      # Local efficiency = sum of inverse distances / possible pairs
      sum(inv_d) / (vcount(subg) * (vcount(subg) - 1))
    })
  }
  
  local_eff <- local_eff_manual(g)
  mean_local_eff <- mean(local_eff, na.rm = TRUE)
  cluster_coeff <- transitivity(g, type = "local", isolates = "zero")
  eig_cent <- eigen_centrality(g, weights = E(g)$weight)$vector
  coreness_values <- coreness(g)
  
  # Compute modularity using community detection
  communities <- cluster_louvain(g, weights = E(g)$weight)
  modularity_score <- modularity(communities)
  
  # Network density
  net_density <- edge_density(g, loops = FALSE)
  
  # Assign community membership to nodes
  V(g)$community <- membership(communities)
  
  # Add centrality values as node attributes
  V(g)$betweenness <- betweenness_centrality
  V(g)$normalized_betweenness <- bet_norm
  V(g)$closeness <- closeness_centrality
  V(g)$degree <- degree_centrality
  V(g)$normalized_degree <- deg_norm
  V(g)$efficiency <- local_eff
  V(g)$eigenvector_centrality <- eig_cent
  V(g)$clustering_coefficient <- cluster_coeff
  V(g)$coreness <- coreness_values
  
  node_list <- data.frame(
    id = V(g)$name,  # Node IDs (this should match the node names in Gephi)
    Label = V(g)$name,
    degree = V(g)$degree,
    normalized_degree = V(g)$normalized_degree,
    betweenness = V(g)$betweenness,
    normalized_betweenness = V(g)$normalized_betweenness,
    closeness = V(g)$closeness,
    community = V(g)$community,
    local_efficiency = V(g)$efficiency,
    eigenvector_centrality = V(g)$eigenvector_centrality,
    local_clustering_coefficient = V(g)$clustering_coefficient,
    Coreness = V(g)$coreness
  )
  
  return(list(
    edges = edge_list, 
    nodes = node_list, 
    linear_relations_df = df_linear,
    modularity = modularity_score,
    mean_local_efficiency = mean_local_eff,
    network_density = net_density))
}

