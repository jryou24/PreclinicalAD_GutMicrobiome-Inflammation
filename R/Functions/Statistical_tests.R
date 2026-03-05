perform_kruskal_dunn <- function(relabund_df, Group = "AD_Status", outdir = getwd(),
                                 output_name){
  
  post_hoc_results <- relabund_df %>%
    group_by(OTU) %>%
    do({
      kruskal_result <- kruskal.test(Relative_abundance ~ .[[Group]], data = .)
      # If the Kruskal-Wallis test is significant (p < 0.05), perform Dunn's post-hoc test
      if(kruskal_result$p.value < 0.05) {
        # Dunn's post-hoc test with Benjamini-Hochberg adjustment
        dunn_result <- dunn.test(.$Relative_abundance, .[[Group]], method = "bh")
        data.frame(
          comparison = dunn_result$comparisons,
          p_value = dunn_result$P.adjusted
        )
      } else {
        data.frame(comparison = NA, p_value = NA)
      }
    }) %>%
    filter(!is.na(comparison)) # Filter out rows with NA comparisons (i.e., non-significant results)
  
  # Prepare data for annotation
  signif_post_hoc <- post_hoc_results %>%
    filter(p_value < 0.1) %>%
    mutate(label = paste(comparison, "\np =", round(p_value, 3))) %>%
    mutate(Group = trimws(str_extract(comparison, "^[^-]+")),
           Group_end = trimws(str_extract(comparison, "[^-]+$")))
  
  print("done")
  write.csv(signif_post_hoc, paste0(outdir,"/",output_name,"_signif_post_hoc.csv"), row.names = F)
  
  return(signif_post_hoc)
 
   
}