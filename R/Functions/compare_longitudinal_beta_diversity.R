within_participant_beta_div_stability <- function(ps_obj, min_interval_in_months = 6, distance_metric = "unifrac",
                                                  distance_metric_string = "Unifrac distance", output_name = "m3_MTG_gt10M_stable_noOther_unif"
){
  # Create pairwise beta diversity matrix
  matrix <- create_pairwise_beta_div_matrix(ps_obj = ps_obj,
                                                      distance_metric = distance_metric,
                                                      distance_metric_string = distance_metric_string,
                                                      output_name = output_name)
  # Filter matrix to within participant comparisons
  matrix_filt <- matrix %>%
    filter(Comparison == "Within participant") 
  
  # Set date columns to date
  matrix_filt$Stool_Date_Collected.x <- as.Date(as.character(matrix_filt$Stool_Date_Collected.x))
  matrix_filt$Stool_Date_Collected.y <- as.Date(as.character(matrix_filt$Stool_Date_Collected.y))
  
  ## Calculate time difference between consecutive samples and only keep comparisons between the baseline and any subsequent samples
  matrix_baseline <- matrix_filt %>% 
    filter(Time_point.x == 1 | Time_point.y == 1) %>%
    group_by(ID.y) %>%
    arrange(Stool_Date_Collected.y) %>%
    mutate(baseline_date = min(Stool_Date_Collected.y, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Time_diff_from_baseline_days = round(abs(as.numeric(difftime(Stool_Date_Collected.x, baseline_date, units = "days"))))) %>%
    mutate(Time_diff_from_baseline_mo = abs(round(time_length(difftime(Stool_Date_Collected.x, baseline_date),"months"),2)))
  
  # Filter out any consecutive samples that are less than 6 months apart from each other
  matrix_baseline_6mo_filt <- matrix_baseline %>%
    arrange(ID.y, Time_diff_from_baseline_mo) %>% 
    group_by(ID.y) %>%
    filter(
      # Step 2: Filter out rows where the difference with the previous row is less than 10
      # Use lag to get the previous Time_diff_from_baseline_mo within each group
      Time_diff_from_baseline_mo - lag(Time_diff_from_baseline_mo, default = -Inf) >= 6
    ) %>%
    ungroup() %>%
    filter(Time_diff_from_baseline_mo >= 6) # Filter out any other samples (tp1 samples) that were collected less than 6 months after the baseline sample
  
  # Creating 1-(beta diversity metric) or stability column
  stability_col_name <- paste0("Stability_", distance_metric)
  matrix_baseline_6mo_filt[[stability_col_name]] <- 1 - matrix_baseline_6mo_filt[[distance_metric]]
  
  return(matrix_baseline_6mo_filt)
}


subset_samples_by_timediff <- function(ps_obj, interval_in_months = 48){
  # Default interval in months will be 48 months (4 years) since our maximum time difference between consecutive samples is around 3 years
  ## Change interval_in_months to desired time interval for beta diversity comparisons between consecutive stool samples 
  library("anytime") 
  library(lubridate)
  library(phyloseq)
  
  # Create output directory if it doesn't exist
  output_path <- paste0(output_path,"_beta_comparisons/")
  if (!file.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  interval_in_months <- as.numeric(interval_in_months)
  
  # obtain metadata table from phyloseq object
  meta <- data.frame(phyloseq::sample_data(ps_obj))

  # changing date columns as dates 
  meta$Date_biomarker_measured <- as.Date(as.character(meta$Date_biomarker_measured)) 
  meta$Stool_Date_Collected <- as.Date(as.character(meta$Stool_Date_Collected)) 
  
  # Extract dataframe of individuals with longitudinal samples and another dataframe of individuals with
  ## single time point samples
  # First calculate the number of entries for each individual
  individual_counts <- meta %>%
    group_by(ID) %>%
    summarise(entry_count = n())

  # Subset longitudinal dataframe
  longitudinal <- meta %>%
    filter(ID %in% individual_counts$ID[individual_counts$entry_count >= 2]) %>%
    mutate(Stool_Date_Collected = as.Date(Stool_Date_Collected, format = "%m/%d/%Y"))
  write.csv(longitudinal, paste0(output_path,deparse(substitute(ps_obj)),"_longitudinal_samples.csv"), 
            row.names = F)
  # Filter out single time point samples and save them in a separate dataframe
  single_time_point <- meta %>%
    filter(ID %in% individual_counts$ID[individual_counts$entry_count == 1]) %>%
    mutate(Stool_Date_Collected = as.Date(Stool_Date_Collected, format = "%m/%d/%Y"))
  write.csv(single_time_point, paste0(output_path,deparse(substitute(ps_obj)),"_single_timepoint_samples.csv"),
            row.names = F)
  
  ## Separate short-term longitudinal from long-term longitudinal 
  # Calculate time differences between consecutive entries
  longitudinal <- longitudinal %>%
    arrange(ID, Stool_Date_Collected) %>%
    group_by(ID) %>%
    mutate(time_diff = as.numeric(difftime(Stool_Date_Collected, lag(Stool_Date_Collected), units = "days")),
           prev_date_of_stool = lag(Stool_Date_Collected)) 
  
  write.csv(longitudinal, paste0(output_path,deparse(substitute(ps_obj)),"_longitudinal_samples_time_diff.csv"),
            row.names = F)
  
  
  # Remove rows with "NA" time differences. These are the first time point samples. 
  longitudinal_pairs <- longitudinal %>%
    filter(!is.na(time_diff)) %>%   # Filter first time point entries since those will have an "NA"
    mutate(time_diff_months = time_diff / 30.4375) # approximate months

  write.csv(longitudinal_pairs, paste0(output_path,deparse(substitute(ps_obj)),"_longitudinal_samples_time_diff_filt.csv"),
            row.names = F)
  
  # Identify rows to keep for short_term
  short_term_ids <- longitudinal_pairs %>%
    filter(time_diff_months <= interval_in_months) %>%
    select(ID, Stool_Date_Collected) %>%
    bind_rows(longitudinal_pairs %>%
                filter(time_diff_months <= interval_in_months) %>%
                select(ID, Stool_Date_Collected = prev_date_of_stool)) %>%
    distinct()
  
  write.csv(short_term_ids, paste0(output_path,deparse(substitute(ps_obj)),"_",
                                   interval_in_months,"_longitudinal_short-term_samples.csv"),
            row.names = F)
  
  # Create short_term dataframe including both entries
  short_term <- longitudinal %>%
    semi_join(short_term_ids, by = c("ID", "Stool_Date_Collected"))
  
  write.csv(short_term, paste0(output_path,deparse(substitute(ps_obj)),"_",
                               interval_in_months,"months_longitudinal_short_term_final.csv"),
            row.names = F)
  # Create long_term dataframe including all entries that are not in short_term
  long_term <- longitudinal %>%
    anti_join(short_term, by = c("ID", "Stool_Date_Collected"))
  
  write.csv(long_term, paste0(output_path,deparse(substitute(ps_obj)),"_",
                              interval_in_months,"months_longitudinal_long_term_final.csv"),
            row.names = F)
    
}

create_pairwise_beta_div_matrix <- function(ps_obj, distance_metric, save_matrix = FALSE, distance_metric_string, output_name){
  # Create output directory if it doesn't exist
  output_path <- paste0(output_path,"_beta_comparisons/")
  if (!file.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  meta <- data.frame(phyloseq::sample_data(ps_obj))
  # Create distance matrix table
  ps.dist <- phyloseq::distance(ps_obj, method=distance_metric)
  ###Looking at pairwise  distances
  ps.dist.pairwise <- as.matrix(ps.dist)
  ps.dist.pairwise2 <- melt(ps.dist.pairwise)
  colnames(ps.dist.pairwise2)[1] <- "Sample_ID1"
  colnames(ps.dist.pairwise2)[2] <- "Sample_ID2"
  colnames(ps.dist.pairwise2)[3] <- distance_metric
  
  # Removing duplicated rows
  ps.dist.pairwise_uniq <- ps.dist.pairwise2[!duplicated(apply(ps.dist.pairwise2, 1, function(x) paste(sort(x), collapse = ""))),]
  # Removing self comparisons with values of zero
  ps.dist.pairwise_uniq <- ps.dist.pairwise_uniq[which(ps.dist.pairwise_uniq$Sample_ID1 != ps.dist.pairwise_uniq$Sample_ID2),] 
  # Merge final distance matrix with metadata
  ps.dist.pairwise_uniq <- merge(ps.dist.pairwise_uniq, meta, by.x = "Sample_ID1", by.y = "Sample_ID")
  ps.dist.pairwise_uniq <- merge(ps.dist.pairwise_uniq, meta, by.x = "Sample_ID2", by.y = "Sample_ID")
  
  ### Comparing within and between participant distances 
  # Adding comparison type column 
  ps.dist.pairwise_uniq$Comparison[ps.dist.pairwise_uniq$ID.x == ps.dist.pairwise_uniq$ID.y] <- "Within participant"
  ps.dist.pairwise_uniq$Comparison[ps.dist.pairwise_uniq$ID.x != ps.dist.pairwise_uniq$ID.y] <- "Between participant"
  
  if(save_matrix == TRUE) {
    write.csv(ps.dist.pairwise_uniq, paste0(output_path, deparse(substitute(ps_obj)), "_", distance_metric_string, "_",
                                            output_name, ".csv"), row.names = F)
  }
  
  ## Prepare for Wilcoxon ranksum tests between AD Status and participant pair groups
  ps.dist.pairwise_uniq$Comparison <- factor(ps.dist.pairwise_uniq$Comparison, levels = c("Between participant","Within participant"))
  
  return(ps.dist.pairwise_uniq)
  
}


subset_phyloseq_obj <- function(ps_obj, target_df){
  ps_obj_meta <- data.frame(phyloseq::sample_data(ps_obj))
  target_sample_list <- as.vector(target_df$Sample_ID)
  target_ps_obj <- phyloseq::subset_samples(ps_obj, phyloseq::sample_names(ps_obj) %in% target_sample_list)
  
  return(target_ps_obj)
}