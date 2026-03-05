### Multiple ecocsystem comparisons assume the samples 
#are the same ie sampling from two different body sites 
# or multiple time points for the same sample
# so to by pass this we will compare correlations across ecosystems
# basically looking for correlated bugs in preclinical vs healthy
# the goal is to see if there are correlations that are status specific letting
# us see if there are ecological shifts that are status specific 
get_upper_tri = function(cormat){
  cormat[lower.tri(cormat)] = NA
  diag(cormat) = NA
  return(cormat)
}
secom_linear_wrapper <- function(ps_obj, seed=42, values=c("healthy"), 
                           output=output_path, method="spearman"){
  set.seed(seed) 
  temp_ps <- ps_obj %>% microViz::ps_filter(AD_Status %in% values)
  
  res = secom_linear(data = list(temp_ps), 
                            assay_name = "counts", tax_level = "Species", 
                            pseudo = 0,lib_cut=0, method =method, 
                            soft = FALSE, thresh_len = 100, n_cv = 10, 
                            thresh_hard = 0, max_p = 0.005)
  corr = res$corr_fl
  cooccur = res$mat_cooccur
  overlap = 10
  corr[cooccur < overlap] = 0
  df_linear = data.frame(get_upper_tri(corr)) %>%
    rownames_to_column("var1") %>%
    pivot_longer(cols = -var1, names_to = "var2", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(value = round(value, 2))
  
  return(df_linear)
}


secom_distance_wrapper <- function(ps_obj, seed=42, values=c("healthy"), 
                                  output=output_path){
  set.seed(seed) 
  temp_ps <- ps_obj %>% microViz::ps_filter(AD_Status %in% values)
  
  res = secom_dist(data = list(ps_obj), 
                   assay_name = "counts", tax_level = "Species", 
                   pseudo = 0,lib_cut=0,
                   thresh_hard = 0, max_p = 0.005)
  corr_dist = res$dcorr_fl
  cooccur_dist = res$mat_cooccur
  overlap = 10
  corr_dist[cooccur_dist < overlap] = 0
  
  df_dist = data.frame(get_upper_tri(corr_dist)) %>%
    rownames_to_column("var1") %>%
    pivot_longer(cols = -var1, names_to = "var2", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(var2 = gsub("\\...", " - ", var2),
           value = round(value, 2))
  
  return(df_dist)
}


secom_comparison <- function(ps_obj, seed=42, value1=c("healthy"), value2=c("preclinical"), 
                output=output_path, method="spearman", lower_cutoff=.5, upper_cutoff=2)
  {
  # cutoffs are the ratio of coefficients
  # the idea is to exclude values greater than the lower cutoff and less than the higher cutoff
  # so a coefficient of .8/.9 would be excluded but -.8/.9 wouldnt be
  # theres a lot of room to play with these, I just included to minimize noise from "preserved"
  # coefficients between two ecosystems being compared
  
  temp1_linear <- secom_linear_wrapper(ps_obj, values=value1) 
  temp2_linear <- secom_linear_wrapper(ps_obj, values=value2)
  
  temp1_dist <- secom_distance_wrapper(ps_obj, values=value1) 
  temp2_dist <- secom_distance_wrapper(ps_obj, values=value2)
  
  merged_linear <- merge(temp1_linear, temp2_linear, by=c("var1", "var2"), all = TRUE) %>% filter(!(value.x == 0 & value.y == 0)) %>% filter(!(value.x/value.y > lower_cutoff & value.x/value.y < upper_cutoff))
  #skipping filtering by ratios since dist coeficcients aren't negative 
  merged_dist <-  merge(temp1_dist, temp2_dist, by=c("var1", "var2"), all = TRUE) %>% filter(!(value.x == 0 & value.y == 0)) 
  
  concatenated_linear <- c(merged_linear$var1, merged_linear$var2)
  unique_values_linear <- unique(concatenated_linear)
  
  concatenated_dist <- c(merged_dist$var1, merged_dist$var2)
  unique_values_dist <- unique(concatenated_dist)
  
  write.csv(unique_values_linear, paste0(output,value1, "_",value2, "_secom_linear_nodes.csv"), row.names = FALSE)
  write.csv(merged_linear, paste0(output,value1, "_",value2, "_secom_Linear_table.csv"), row.names = FALSE)
  
  write.csv(unique_values_dist, paste0(output,value1, "_",value2, "_secom_dist_nodes.csv"), row.names = FALSE)
  write.csv(merged_dist, paste0(output,value1, "_",value2, "_secom_dist_table.csv"), row.names = FALSE)
  
}




