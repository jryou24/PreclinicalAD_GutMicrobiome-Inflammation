subset_sample_closest_to_biomarker_measure <- function(ps_obj) {
  library("anytime") 
  library(lubridate)
  library(phyloseq)
  
  # obtain metadata table from phyloseq object
  meta <- data.frame(phyloseq::sample_data(ps_obj))
  
  # Filter metadata table by shortest time difference between stool date and biomarker measurement date
  final_df <- meta %>%
    group_by(ID) %>%
    filter(Days_since_biomarker_assessment == min(Days_since_biomarker_assessment)) %>%
    ungroup() 

  return(final_df)

}