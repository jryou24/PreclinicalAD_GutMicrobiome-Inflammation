

perform_mtx_model_DE_all <- function(mtx_obj, mtg_obj, model, fixed_vars, reference, random_vars, output_name, 
                                     rna_dna_filt = "semistrict",  core=8, filter_vector = c("Healthy", "Preclinical", "Symptomatic")){
  
  output_folder=paste0(output_path,"_MTXmodel/",deparse(substitute(mtx_obj)),"_",output_name, "_", rna_dna_filt)
  dir.create(output_folder, recursive = TRUE)
  # To fit a negative binomial model, count data is required. Here we transform
  # relative abundances to integer count data, preserving accuracy to 0.00001%.  
  mtg_obj <- mtg_obj %>% microViz::ps_filter(AD_Status %in% filter_vector)
  mtx_obj <- mtx_obj %>% microViz::ps_filter(AD_Status %in% filter_vector)
  
  mtx_obj_int <- transform_sample_counts(mtx_obj, function(x) trunc(x*100000))
  mtg_obj_int <- transform_sample_counts(mtg_obj, function(x) trunc(x*100000))
  
  # Access metadata and OTU tables 
  meta.de <- data.frame(sample_data(mtx_obj))
  mtx_input.int <- data.frame(otu_table(mtx_obj_int), check.names = FALSE) 
  mtg_input.int <- data.frame(otu_table(mtg_obj_int), check.names = FALSE) 
  
  print("Reading of files and data transformations complete")
  print(dim(mtx_input.int))
  print(dim(mtg_input.int))
  print("Starting DE analysis using MTXmodel")
  
  fixed.vars <- as.vector(fixed_vars)
  
  fit_data <- MTXmodel(
    mtx_input.int, meta.de, output_folder, transform = "NONE",
    analysis_method= model, 
    fixed_effects = fixed.vars,
    random_effects = random_vars,
    reference = reference,
    normalization = 'NONE',
    standardize = FALSE,
    input_dnadata = mtg_input.int,
    rna_dna_flt = rna_dna_filt,
    cores = core)
  
}


maaslin_wrapper <- function(ps_obj, reference_vector, random_effects_vector,
                            fixed_effects_vector=c("sex", "age","DIABETES", "BMI", "cancer", "stool_date_cont", "AD_Status"), 
                            output=output_path, filter_vector = c("healthy", "preclinical", "symptomatic"), norm_method = "NONE", anal_method='NEGBIN', 
                            transf_method = "NONE", core=8){
  dir.create(paste0(output_path,"_maaslin/",deparse(substitute(ps_obj))), recursive = TRUE)
  # Docs and code is from Aura 
  # Use taxa abundances from phyloseq object ps2.filt, which has already been 
  # filtered to omit taxa with < 0.1% mean relative abundance across all samples. 
  # To fit a negative binomial model, count data is required. Here we transform
  # relative abundances to integer count data, preserving accuracy to 0.00001%.  
  #temp_ps <- ps_obj %>% microViz::ps_filter(AD_Status %in% filter_vector)
  temp_ps <- transform_sample_counts(ps_obj, function(x) trunc(x*100000))
  
  species.ps2<- data.frame(otu_table(temp_ps), check.names = FALSE) 
  metadata.ps2 <- data.frame(sample_data(temp_ps), check.names = FALSE)
  
  
  
  #A run Maaslin at species level, with pre-clinical status --
  fit_Maaslin2 <- Maaslin2(
    input_data = species.ps2,
    input_metadata = metadata.ps2,
    normalization = norm_method,      # Already normalized (count-transformed relabund)
    standardize = 'TRUE',        # We'd like to z-score continuous metadata
    min_prevalence = 0,          # No filtering on prevalence
    min_abundance = 0,           # Data already filtered on abundance
    transform = transf_method,          # Further transform not needed for NEGBIN model
    analysis_method = anal_method,
    max_significance = 0.05,     # For q-values (BH-adjusted p-values). Default is 0.25.
    output = paste0(output_path,"_maaslin/",deparse(substitute(ps_obj))),
    fixed_effects = fixed_effects_vector, 
    random_effects = random_effects_vector,
    reference=reference_vector,
    cores = core, 
    plot_heatmap = FALSE,
    plot_scatter = FALSE

  )
}


maaslin_linear_wrapper <- function(ps_obj, reference_vector, random_effects_vector, prev=0, abun=0,
                            fixed_effects_vector=c("sex", "age","DIABETES", "BMI", "cancer", "stool_date_cont", "AD_Status"), 
                            output=output_path, filter_vector = c("healthy", "preclinical", "symptomatic"), anal_method='LM',correction="BH",  core=8,norm='NONE', transform='NONE'){
  dir.create(paste0(output_path,"_maaslin_linear/",deparse(substitute(ps_obj))), recursive = TRUE)
  # Docs and code is from Aura 
  # Use taxa abundances from phyloseq object ps2.filt, which has already been 
  # filtered to omit taxa with < 0.1% mean relative abundance across all samples. 
  # To fit a negative binomial model, count data is required. Here we transform
  # relative abundances to integer count data, preserving accuracy to 0.00001%.  
  #temp_ps <- ps_obj %>% microViz::ps_filter(AD_Status %in% filter_vector)
  #temp_ps <- transform_sample_counts(temp_ps, function(x) trunc(x*100000))
  
  species.ps2<- data.frame(otu_table(ps_obj), check.names = FALSE) 
  metadata.ps2 <- data.frame(sample_data(ps_obj), check.names = FALSE)
  
  
  
  #A run Maaslin at species level, with pre-clinical status --
  fit_Maaslin2 <- Maaslin2(
    input_data = species.ps2,
    input_metadata = metadata.ps2,
    normalization = norm,      # Already normalized (count-transformed relabund)
    standardize = 'TRUE',        # We'd like to z-score continuous metadata
    min_prevalence = prev,          # No filtering on prevalence
    min_abundance = abun,           # Data already filtered on abundance
    transform = transform,          # Further transform not needed for NEGBIN model
    analysis_method = anal_method,
    max_significance = 0.05,     # For q-values (BH-adjusted p-values). Default is 0.25.
    output = paste0(output_path,"_maaslin_linear/",deparse(substitute(ps_obj))),
    fixed_effects = fixed_effects_vector, 
    random_effects = random_effects_vector,
    reference=reference_vector,
    cores = core,
    correction=correction
  )
}

