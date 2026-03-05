
make_phyloseq_object <- function(path_to_metaphlan_txt, path_to_metadata, tree_object, skip=1, cat_vars=c('ID', 'Sample_ID', 'Time_point', 'AD_Status', 'sex', 'race', 'cdr', 'apoe', 'CBSTROKE','DIABETES','HYPERTEN','HYPERCHO','CVHATT','THYROID', 'DEP2YRS','ALCOHOL','INDEPEND','RESIDENC','cancer','APP','APPSEV','Tracer_used_for_classification'), num_vars, abund_threshold, columns_to_drop=c(), compute_mean=FALSE, sum_cutoff=0.01, subset_samples=c()){
  
  # Read in m4 table
  ## Metaphlan files will often have a first line with metadata, skip=1 should let you read in the file without needing to delete the line

  mphlan <- data.frame(read.table(path_to_metaphlan_txt, header = T, sep="\t", skip=skip))
  print("Done reading metaphlan table")
  
  if (typeof(path_to_metadata) == "character"){
    print("Metadata path provided, reading it in")
  metadata <- data.frame(read_excel(path_to_metadata))
  # Convert character coded "NA" to actual NA
  metadata <- metadata %>%
    mutate(across(where(is.character), ~ na_if(., "NA")))
  }
  else{
    metadata <- data.frame(path_to_metadata)
  }

  print("Done reading metadata table")
  
  # Filter out samples that had 100% unclassified reads (probably due to low sequencing depth)
  #unclassified_row <- which(mphlan$clade_name == "UNCLASSIFIED")
  #mphlan <- mphlan[, mphlan[unclassified_row, ] != 100]
  
  ### Processing for phyloseq takes 3 steps 1) species data 2) taxa table 3) metadata
  #1) Species Data 
  rownames(mphlan)<-mphlan$clade_name
  mphlan$clade_name<- NULL
  # keep rows with species info
  mphlan.species <- mphlan[grepl('s__', rownames(mphlan)), ]
  # Drop rows with SGB (species level genome bins)
  mphlan.species <- mphlan.species[!grepl('t__', rownames(mphlan.species)), ]
  mphlan.species <- mphlan.species[, !(names(mphlan.species) %in% columns_to_drop)]
  
  if(compute_mean==TRUE){
    
    colnames_species <- colnames(mphlan.species)
    # strip only "AD_id#" from "AD_id#_tp" from colnames
    colnames_species_no_tp <- gsub("(.+?_.+?)_.*","\\1", colnames_species)
    
    unique_ind <- unique(colnames_species_no_tp)
    unique_ind <- as.vector(unique_ind)
    
    #
    temp_species_df <- mphlan.species
    colnames(temp_species_df) <- colnames_species_no_tp
    mean_species <- data.frame(rownames(mphlan.species))
    colnames(mean_species) = "Species"
    
    for (ind in unique_ind) {
      
      df <- data.frame(temp_species_df) %>% 
        dplyr::select(starts_with(ind))
      df$mean = apply(df, 1, mean, na.rm = T) # calculate mean across rows
      mean_species = cbind(mean_species, mean = df$mean)
      names(mean_species)[names(mean_species) == 'mean'] <- ind
    }
    rownames(mean_species) <- rownames(mphlan.species)
    colnames(mean_species) <- paste(colnames(mean_species),"01",sep="_") 
    mean_species$Species_01<-NULL
    mphlan.species <- mean_species
    mphlan.species <- as.matrix(mphlan.species)
    
    
  }
  else{
    mphlan.species <- as.matrix(mphlan.species)
    
  }
  print(dim(mphlan.species))
  #2) Taxa table 
  # Create phyloseq taxa table
  print("Making taxa table")
  
  species <- data.frame(Names = rownames(mphlan.species))
  species <- data.frame(do.call('rbind', strsplit(as.character(species$Names),'|', fixed=TRUE)))
  rownames(species) <- rownames(mphlan.species )
  colnames(species) <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
  print(dim(species))
  species <- as.matrix(species)
  #3) Metadata 
  # Filter out samples not found in the metaphlan results and adjust factors
  rows_metadata <- rownames(metadata)[(metadata$Sample_ID %in% colnames(mphlan.species))]
  metadata <- metadata[rownames(metadata) %in% rows_metadata, ]
  metadata <- data.frame(metadata)
  rownames(metadata) <- metadata$Sample_ID 
  # Set categorical variables to factor
  metadata[,cat_vars] <- lapply(metadata[,cat_vars], factor)
  # Set numerical variables to numeric
  metadata[,num_vars] <- lapply(metadata[,num_vars], as.numeric)
  print(dim(metadata))
  
  print("Differences between metadata rows and metaphlan column names")
  print(setdiff(colnames(mphlan.species), rownames(metadata)))
  ps <- phyloseq(otu_table(mphlan.species, taxa_are_rows = TRUE),
                 sample_data(metadata),
                 tax_table(species), 
                 phy_tree(tree_object)) 
  ps.filt <- prune_taxa(taxa_sums(ps)>sum_cutoff, ps)
  ps.filt <- filter_taxa(ps.filt, function(x) mean(x) > abund_threshold, TRUE)
  ps.filt <- transform_sample_counts(ps.filt, function(x) (x / sum(x))*100)
  
  return(ps.filt)
}


make_phyloseq_pathway_object <- function(path_to_featuretable_txt, path_to_metadata, cat_vars, num_vars,abundance_threshold, prevalence_threshold,columns_to_drop=c(), filter_species_associations=FALSE, filter_by_abundance = FALSE, filter_by_prevalence_abundance = FALSE, filter_by_prevalence = TRUE){
  
  # Read in humann feature table
  humann <- data.frame(read_tsv(path_to_featuretable_txt))
  ### Processing for phyloseq takes 3 steps 1) species data 2) taxa table 3) metadata
  #1) Species Data 
  ## Rename name of first column 
  colnames(humann)[1] <- "Pathway"
  humann$Pathway <- as.character(humann$Pathway)
  
  # If filter_species_associations is set to TRUE, omit rows that break down the species associations (contain '|' delimiter)
  if(filter_species_associations==TRUE){
    humann <- humann %>%
      filter(!str_detect(Pathway, "\\|"))
  }
  else{
    # If not, omit rows that DON'T break down the species associations so we only retain pathways/KOs that have them
    humann <- humann %>%
      filter(str_detect(Pathway, "\\|"))
  }
  # Replace unwanted characters from pathway names with underscores
  #[^a-zA-Z0-9\s] matches any character that is not alphanumeric (a-z, A-Z, 0-9) or whitespace (\s).
  humann$Pathway <- gsub("[^[:alnum:]\\s]","_", humann$Pathway)


  # Filter out unintegrated/unclassified/unmapped features 
  humann <- humann %>%
    filter(!str_detect(Pathway, "UNMAPPED")) %>%
    filter(!str_detect(Pathway, "UNINTEGRATED")) %>%
    filter(!str_detect(Pathway, "unclassified")) %>%
    filter(!str_detect(Pathway, "UNGROUPED"))

  # obtain number of samples
  nsamples = ncol(humann)
  
  # Set pathway column as rownames 
  rownames(humann)<-humann$Pathway
  humann$Pathway<- NULL
  # Check pathway relative abundances sum to 1 (unlike taxa, which sum to 100)
  humann.sums <- colSums(humann) 
  min_abund <- min(humann.sums) 
  print(paste("The minimum pathway abundance is:", min_abund))
  print("Filtering samples that failed humann profiling or sequencing")
  # filter out columns with zero sum
  humann <- humann %>% dplyr::select_if(~sum(.) != 0)

  print(paste("The dimensions of the pathway abundance table are:", dim(humann))) 
  humann <- as.matrix(humann)
  
  #2) Taxa table 
  # Create phyloseq taxa table
  print("Making taxa table")
  pathways <- data.frame(Pathway=rownames(humann))
  rownames(pathways) <- pathways$Pathway
  pathways <- as.matrix(pathways)
  print(dim(pathways))
  
  print(paste("number of samples:",nsamples))
  
  #3) Metadata 
  metadata <- data.frame(read_xlsx(path_to_metadata))
  rownames(metadata) <- metadata$Sample_ID 
  # Set categorical variables to factor
  metadata[,cat_vars] <- lapply(metadata[,cat_vars], factor)
  # Set numerical variables to numeric
  metadata[,num_vars] <- lapply(metadata[,num_vars], as.numeric)
  
  print("Differences between metadata rows and metaphlan column names")
  print(setdiff(colnames(humann), rownames(metadata)))
  
  print("Checking dimensions before phyloseq object creation")
  print(dim(humann))
  print(dim(metadata))
  print(dim(pathways))
  
  # Create phyloseq object
  ps <- phyloseq(otu_table(humann, taxa_are_rows = TRUE),
                 sample_data(metadata),
                 tax_table(pathways)) 
  print("Phyloseq object successfully created!")

  # Renormalize relative abundances 
  ps  <- transform_sample_counts(ps, function(x) (x / sum(x))*100 )
  temp_otu <- data.frame(otu_table(ps), check.names = FALSE) 
  col_sums <- colSums(temp_otu)
  print(col_sums)
  
  # Filter lowly abundant pathways depending on the condition
  if(filter_by_prevalence_abundance==TRUE){
    print("Filtering OTU table by prevalence and abundance")
    ps.filt <- filter_taxa(ps, function(x)(sum(x > abundance_threshold) > (nsamples*prevalence_threshold)), prune = TRUE) 
    
  }
  if(filter_by_abundance==TRUE){
    print("Filtering OTU table by abundance only")
    ps.filt <- filter_taxa(ps, function(x) mean(x) > abundance_threshold, prune = TRUE) 
    
  }
  if(filter_by_prevalence==TRUE){
    print(paste0("Filtering OTU table by prevalence > ", prevalence_threshold))
    ps.filt <- filter_taxa(ps, function(x) (mean(x>0) > prevalence_threshold), prune = TRUE) 
  }
  # Renormalize relative abundances 
  ps.filt  <- transform_sample_counts(ps.filt, function(x) (x / sum(x))*100 )
  
  return(ps.filt)
}


make_phyloseq_object <- function(path_to_metaphlan_txt, path_to_metadata, tree_object, skip=1, cat_vars=c('ID', 'Sample_ID', 'Time_point', 'AD_Status', 'sex', 'race', 'cdr', 'apoe', 'CBSTROKE','DIABETES','HYPERTEN','HYPERCHO','CVHATT','THYROID', 'DEP2YRS','ALCOHOL','INDEPEND','RESIDENC','cancer','APP','APPSEV','Tracer_used_for_classification'), num_vars, abund_threshold, columns_to_drop=c(), compute_mean=FALSE, sum_cutoff=0.01, subset_samples=c()){
  
  # Read in m4 table
  ## Metaphlan files will often have a first line with metadata, skip=1 should let you read in the file without needing to delete the line

  mphlan <- data.frame(read.table(path_to_metaphlan_txt, header = T, sep="\t", skip=skip))
  print("Done reading metaphlan table")
  
  if (typeof(path_to_metadata) == "character"){
    print("Metadata path provided, reading it in")
  metadata <- data.frame(read_excel(path_to_metadata))
  # Convert character coded "NA" to actual NA
  metadata <- metadata %>%
    mutate(across(where(is.character), ~ na_if(., "NA")))
  }
  else{
    metadata <- data.frame(path_to_metadata)
  }

  print("Done reading metadata table")
  
  # Filter out samples that had 100% unclassified reads (probably due to low sequencing depth)
  #unclassified_row <- which(mphlan$clade_name == "UNCLASSIFIED")
  #mphlan <- mphlan[, mphlan[unclassified_row, ] != 100]
  
  ### Processing for phyloseq takes 3 steps 1) species data 2) taxa table 3) metadata
  #1) Species Data 
  rownames(mphlan)<-mphlan$clade_name
  mphlan$clade_name<- NULL
  # keep rows with species info
  mphlan.species <- mphlan[grepl('s__', rownames(mphlan)), ]
  # Drop rows with SGB (species level genome bins)
  mphlan.species <- mphlan.species[!grepl('t__', rownames(mphlan.species)), ]
  mphlan.species <- mphlan.species[, !(names(mphlan.species) %in% columns_to_drop)]
  
  if(compute_mean==TRUE){
    
    colnames_species <- colnames(mphlan.species)
    # strip only "AD_id#" from "AD_id#_tp" from colnames
    colnames_species_no_tp <- gsub("(.+?_.+?)_.*","\\1", colnames_species)
    
    unique_ind <- unique(colnames_species_no_tp)
    unique_ind <- as.vector(unique_ind)
    
    #
    temp_species_df <- mphlan.species
    colnames(temp_species_df) <- colnames_species_no_tp
    mean_species <- data.frame(rownames(mphlan.species))
    colnames(mean_species) = "Species"
    
    for (ind in unique_ind) {
      
      df <- data.frame(temp_species_df) %>% 
        dplyr::select(starts_with(ind))
      df$mean = apply(df, 1, mean, na.rm = T) # calculate mean across rows
      mean_species = cbind(mean_species, mean = df$mean)
      names(mean_species)[names(mean_species) == 'mean'] <- ind
    }
    rownames(mean_species) <- rownames(mphlan.species)
    colnames(mean_species) <- paste(colnames(mean_species),"01",sep="_") 
    mean_species$Species_01<-NULL
    mphlan.species <- mean_species
    mphlan.species <- as.matrix(mphlan.species)
    
    
  }
  else{
    mphlan.species <- as.matrix(mphlan.species)
    
  }
  print(dim(mphlan.species))
  #2) Taxa table 
  # Create phyloseq taxa table
  print("Making taxa table")
  
  species <- data.frame(Names = rownames(mphlan.species))
  species <- data.frame(do.call('rbind', strsplit(as.character(species$Names),'|', fixed=TRUE)))
  rownames(species) <- rownames(mphlan.species )
  colnames(species) <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
  print(dim(species))
  species <- as.matrix(species)
  #3) Metadata 
  # Filter out samples not found in the metaphlan results and adjust factors
  rows_metadata <- rownames(metadata)[(metadata$Sample_ID %in% colnames(mphlan.species))]
  metadata <- metadata[rownames(metadata) %in% rows_metadata, ]
  metadata <- data.frame(metadata)
  rownames(metadata) <- metadata$Sample_ID 
  # Set categorical variables to factor
  metadata[,cat_vars] <- lapply(metadata[,cat_vars], factor)
  # Set numerical variables to numeric
  metadata[,num_vars] <- lapply(metadata[,num_vars], as.numeric)
  print(dim(metadata))
  
  print("Differences between metadata rows and metaphlan column names")
  print(setdiff(colnames(mphlan.species), rownames(metadata)))
  ps <- phyloseq(otu_table(mphlan.species, taxa_are_rows = TRUE),
                 sample_data(metadata),
                 tax_table(species), 
                 phy_tree(tree_object)) 
  ps.filt <- prune_taxa(taxa_sums(ps)>sum_cutoff, ps)
  ps.filt <- filter_taxa(ps.filt, function(x) mean(x) > abund_threshold, TRUE)
  ps.filt <- transform_sample_counts(ps.filt, function(x) (x / sum(x))*100)
  
  return(ps.filt)
}


make_phyloseq_pathway_object <- function(path_to_featuretable_txt, path_to_metadata, cat_vars, num_vars,abundance_threshold, prevalence_threshold,columns_to_drop=c(), filter_species_associations=FALSE, filter_by_abundance = FALSE, filter_by_prevalence_abundance = FALSE, filter_by_prevalence = TRUE){
  
  # Read in humann feature table
  humann <- data.frame(read_tsv(path_to_featuretable_txt))
  ### Processing for phyloseq takes 3 steps 1) species data 2) taxa table 3) metadata
  #1) Species Data 
  ## Rename name of first column 
  colnames(humann)[1] <- "Pathway"
  humann$Pathway <- as.character(humann$Pathway)
  
  # If filter_species_associations is set to TRUE, omit rows that break down the species associations (contain '|' delimiter)
  if(filter_species_associations==TRUE){
    humann <- humann %>%
      filter(!str_detect(Pathway, "\\|"))
  }
  else{
    # If not, omit rows that DON'T break down the species associations so we only retain pathways/KOs that have them
    humann <- humann %>%
      filter(str_detect(Pathway, "\\|"))
  }
  # Replace unwanted characters from pathway names with underscores
  #[^a-zA-Z0-9\s] matches any character that is not alphanumeric (a-z, A-Z, 0-9) or whitespace (\s).
  humann$Pathway <- gsub("[^[:alnum:]\\s]","_", humann$Pathway)


  # Filter out unintegrated/unclassified/unmapped features 
  humann <- humann %>%
    filter(!str_detect(Pathway, "UNMAPPED")) %>%
    filter(!str_detect(Pathway, "UNINTEGRATED")) %>%
    filter(!str_detect(Pathway, "unclassified")) %>%
    filter(!str_detect(Pathway, "UNGROUPED"))

  # obtain number of samples
  nsamples = ncol(humann)
  
  # Set pathway column as rownames 
  rownames(humann)<-humann$Pathway
  humann$Pathway<- NULL
  # Check pathway relative abundances sum to 1 (unlike taxa, which sum to 100)
  humann.sums <- colSums(humann) 
  min_abund <- min(humann.sums) 
  print(paste("The minimum pathway abundance is:", min_abund))
  print("Filtering samples that failed humann profiling or sequencing")
  # filter out columns with zero sum
  humann <- humann %>% dplyr::select_if(~sum(.) != 0)

  print(paste("The dimensions of the pathway abundance table are:", dim(humann))) 
  humann <- as.matrix(humann)
  
  #2) Taxa table 
  # Create phyloseq taxa table
  print("Making taxa table")
  pathways <- data.frame(Pathway=rownames(humann))
  rownames(pathways) <- pathways$Pathway
  pathways <- as.matrix(pathways)
  print(dim(pathways))
  
  print(paste("number of samples:",nsamples))
  
  #3) Metadata 
  metadata <- data.frame(read_xlsx(path_to_metadata))
  rownames(metadata) <- metadata$Sample_ID 
  # Set categorical variables to factor
  metadata[,cat_vars] <- lapply(metadata[,cat_vars], factor)
  # Set numerical variables to numeric
  metadata[,num_vars] <- lapply(metadata[,num_vars], as.numeric)
  
  print("Differences between metadata rows and metaphlan column names")
  print(setdiff(colnames(humann), rownames(metadata)))
  
  print("Checking dimensions before phyloseq object creation")
  print(dim(humann))
  print(dim(metadata))
  print(dim(pathways))
  
  # Create phyloseq object
  ps <- phyloseq(otu_table(humann, taxa_are_rows = TRUE),
                 sample_data(metadata),
                 tax_table(pathways)) 
  print("Phyloseq object successfully created!")

  # Renormalize relative abundances 
  ps  <- transform_sample_counts(ps, function(x) (x / sum(x))*100 )
  temp_otu <- data.frame(otu_table(ps), check.names = FALSE) 
  col_sums <- colSums(temp_otu)
  print(col_sums)
  
  # Filter lowly abundant pathways depending on the condition
  if(filter_by_prevalence_abundance==TRUE){
    print("Filtering OTU table by prevalence and abundance")
    ps.filt <- filter_taxa(ps, function(x)(sum(x > abundance_threshold) > (nsamples*prevalence_threshold)), prune = TRUE) 
    
  }
  if(filter_by_abundance==TRUE){
    print("Filtering OTU table by abundance only")
    ps.filt <- filter_taxa(ps, function(x) mean(x) > abundance_threshold, prune = TRUE) 
    
  }
  if(filter_by_prevalence==TRUE){
    print(paste0("Filtering OTU table by prevalence > ", prevalence_threshold))
    ps.filt <- filter_taxa(ps, function(x) (mean(x>0) > prevalence_threshold), prune = TRUE) 
  }
  # Renormalize relative abundances 
  ps.filt  <- transform_sample_counts(ps.filt, function(x) (x / sum(x))*100 )
  
  return(ps.filt)
}





