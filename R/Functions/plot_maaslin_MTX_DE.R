plot_coefficient_prevalence_maaslin <- function(maaslin_output, ps_obj, num_samples=188, var, var_levels, values=c("healthy", "preclinical", "symptomatic"),
                                        reference, output_name, colors_dict){
  # Create output directory if it doesn't exist
  output_path <- paste0(output_path,"_maaslin/", deparse(substitute(ps_obj)))
  if (!file.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  sample_cutoff <- num_samples*.10
  #temp_ps <- ps_obj %>% microViz::ps_filter(var %in% values)
  
  maaslin_df <- read_tsv(maaslin_output) %>% filter(metadata==var) %>% mutate(abs_coef=abs(coef)) %>% mutate(qval = p.adjust(pval, method = 'BH')) %>% filter(N.not.0>sample_cutoff) %>% mutate(feature=gsub(".*s__", "", feature)) %>% mutate(value=as.factor(value)) 
  
  median <- median(maaslin_df$abs_coef)
  print(median)
  
  write.csv(maaslin_df, paste0(output_path,"_",var,"_reference_",reference,"_",output_name,"_significant_results.csv"), row.names = F)
  
  print(paste("variable of interest is ", var))
  print("Starting loop of variable levels")
  # Subset to features significantly associated with var, detected in 
  # at least 15% of samples, and with magnitude of the coefficient > median
  var_levels <- as.vector(var_levels)
  print(var_levels)
  # Create empty list to store coefficient and prevalence plots
  return_list = list()
  for (level in var_levels) {
    print(paste("Var is ", var))
    print(paste("Processing ", level))
    name_string <- paste0("diff_feature_",var,"_",level,"_vs_",reference,"_", output_name)
    
    filename <- name_string
    filename <- subset(maaslin_df, value == level  & N.not.0 > sample_cutoff & abs(coef) > median)
    filename <- filename %>% arrange(desc(coef))
    
    #obtain min and max coefficient values (needed for plotting later)
    min <- min(filename$coef)
    print(min)
    max <- max(filename$coef)
    print(max)
    
    ### Generate coefficient plot ####
    # Calculate 95% CI for coefficients
    filename <- filename %>% mutate(pathways = fct_reorder(feature, coef))
    
    filename$coefmax <- filename$coef + filename$stderr
    filename$coefmin <- filename$coef - filename$stderr
    
    filename$CI95upper <- filename$coef + 1.96*filename$stderr
    filename$CI95lower <- filename$coef - 1.96*filename$stderr
    
    # Plot coefficients
    plotname <- paste0(name_string,"_coeff_plot.png")
    x_min_lim <- (floor(min / 5) * 6) 
    print(x_min_lim)
    x_max_lim <- (ceiling(max / 5) * 6) 
    
    write.csv(filename, paste0(output_path,"_",var,"_reference_",reference,"_",output_name,"_significant_results_subset_for_coef_plot.csv"), row.names = F)
    
    print("Plotting coefficient plot")
    
    # create list of colors from
    print(reference)
    print(level)
    color1 <- colors_dict[[reference]]
    color2 <- colors_dict[[level]]
    colors_list <- c(color2, color1)
    print(colors_list)
    
    
    coeff_p <- ggplot(filename, aes(x=coef, y=pathways,color=as.factor(coef<0)))+
      geom_point(size = 3)+
      geom_errorbarh(aes(xmax = CI95upper, xmin = CI95lower, height =0.1))+
      geom_vline(xintercept = 0, linetype='dotted', color='grey', size=1.5)+
      scale_color_manual(values=colors_list)+
      theme_linedraw()+
      theme(axis.text.y = element_text(size=18),
            axis.text.x = element_text(size=18),
            axis.title.x = element_text(size=17),
            axis.title.y = element_blank(),
            legend.text = element_text(size = 16),
            legend.position = "bottom",
            plot.title = element_text(size=15),
            panel.border = element_rect(fill=NA, colour="black", size=3))+
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      xlab('Coefficient')
    #save_plot(paste0(output_path,"_",plotname), coeff_p, base_height=12, base_width=11)
    coeff_p
    ggsave(paste0(output_path,"_",plotname), plot = coeff_p, bg = "transparent", width = 11, height = 12)
    
    #Add coefficient plot to return list
    return_list[[paste0("Coeff_Plot_", level)]] <- coeff_p
    print("Done")
    
    #### Generate prevalence plot ###
    ### Access absolute pathway abundance data 
    abs_abund <- psmelt(ps_obj)
    write.csv(abs_abund, paste0(output_path,"_",name_string,"_abs_abund_original.csv"), 
              row.names = F)
    abs_abund <- abs_abund %>% mutate(OTU=gsub(".*s__", "", OTU))
    diff_feat_abs_abund <- subset(abs_abund, OTU %in% filename$feature)
    write.csv(diff_feat_abs_abund, paste0(output_path,"_",name_string,"_abs_abund_subset_for_prev_plot.csv"), 
              row.names = F)
    # Filter absolute pathway abundance data by variable and level of interest
    list <- c(reference, level)
    print(list)
    diff_feat_abs_abund <- diff_feat_abs_abund %>%
      filter(!!sym(var) %in% list)
    
    # Re-order feature abs abundances according to coefficient value
    # This order is preserved in filename, which was previously ordered: 
    order <- data.frame(OTU = unique(filename$feature))
    # sort = F makes sure that the order is based on order_index df which is ordered by coefficient value
    abs_abund.df <- merge(order, diff_feat_abs_abund, by="OTU", sort = F)
    abs_abund.df$OTU <- as.factor(abs_abund.df$OTU)
    abs_abund.df <- abs_abund.df %>% mutate(OTU = fct_relevel(OTU, filename$feature))
    
    # For each significant feature in each sample, set AbunNotZero == 1 if its 
    # abundance is greater than 0. 
    abs_abund.df$AbunNotZero <- 1*(abs_abund.df$Abundance > 0)
    
    # Summarize pathway prevalence by variable of interest
    abs_abund.df.n0 <- abs_abund.df %>% 
      group_by(OTU, !!sym(var)) %>%
      summarise(N = n(), SumNotZero = sum(AbunNotZero)) %>%
      mutate(PercentNotZero = 100*(SumNotZero / N)) %>%
      ungroup()
    
    
    write.csv(abs_abund.df.n0, paste0(output_path,"_",name_string,"_prevalence.csv"), 
              row.names = F)
    
    # Plot prevalence by variable of interest
    print("Plotting prevalence plot")
    
    plotname2 <- paste0(name_string,"_prev_plot.png")
    
    prev_p <- ggplot(abs_abund.df.n0, aes(x = !!sym(var),y=OTU, fill=PercentNotZero))+
      geom_tile()+
      scale_fill_distiller(palette ='RdBu', type='seq', breaks = c(0,50,100), 'Prevalence (%)')+
      theme_classic()+
      theme(axis.text.y = element_text(size=18),
            axis.text.x = element_text(size=18, angle=45, hjust=1),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size=17),
            legend.text = element_text(size = 16),
            panel.border = element_rect(fill=NA, colour="black", size=3))+
      scale_y_discrete(limits=rev)
    #save_plot(paste0(output_path,"_",plotname2), prev_p, base_height=12, base_width=11)
    ggsave(paste0(output_path,"_",plotname2), plot = prev_p, bg = "transparent", width = 11, height = 12)
    
    # Add prevalence plot to the return list
    return_list[[paste0("Prev_Plot_", level)]] <- prev_p
    print("Done")
    
  }
  return(return_list)
}


plot_coefficient_prevalence <- function(sign_features_df, ps_obj, var, var_levels, 
                                        reference, output_name, colors_dict){
  
  output_folder=paste0(output_path,"_MTXmodel_plots/",deparse(substitute(ps_obj)),"_",output_name,"/")
  dir.create(output_folder, recursive = TRUE)
  
  #Filter significant results table to obtain pathways associated to variable of interest
  results <- sign_features_df %>% filter(metadata == var)
  #re-adjust p-value after filtering using BH method
  results$qval <- p.adjust(results$pval, method = 'BH')
  #create an absolute coefficient column to determine median coefficient value
  results$abs_coef <- abs(results$coef)
  #obtain median
  median <- median(results$abs_coef)
  print(median)
  
  write.csv(results, paste0(output_folder,var,"_significant_results.csv"), row.names = F)
  
  print(paste("variable of interest is ", var))
  print("Starting loop of variable levels")
  # Subset to features significantly associated with var, detected in 
  # at least 15% of samples, and with magnitude of the coefficient > median
  var_levels <- as.vector(var_levels)
  for (level in var_levels) {
    print(paste("Var is ", var))
    print(paste("Processing ", level))
    name_string <- paste0("diff_feature_",var,"_",level)
    
    filename <- name_string
    filename <- subset(results, value == level  & N.not.0 > 53 & abs(coef) > median)
    filename <- filename %>% arrange(desc(coef))
    
    #obtain min and max coefficient values (needed for plotting later)
    min <- min(filename$coef)
    max <- max(filename$coef)
    
    ### Generate coefficient plot ####
    # Calculate 95% CI for coefficients
    filename <- filename %>% mutate(pathways = fct_reorder(feature, coef))
    
    filename$coefmax <- filename$coef + filename$stderr
    filename$coefmin <- filename$coef - filename$stderr
    
    filename$CI95upper <- filename$coef + 1.96*filename$stderr
    filename$CI95lower <- filename$coef - 1.96*filename$stderr
    
    # Plot coefficients
    plotname <- paste0(name_string,"_coeff_plot.png")
    x_min_lim <- (ceiling(min / 10) * 10) - 20
    x_max_lim <- (ceiling(max / 10) * 10) + 20
    
    print("Plotting coefficient plot")
    
    # create list of colors from
    print(reference)
    print(level)
    color1 <- colors_dict[[reference]]
    color2 <- colors_dict[[level]]
    colors_list <- c(color2, color1)
    print(colors_list)
    
    coeff_p <- ggplot(filename, aes(x=coef, y=pathways,color=as.factor(coef<0)))+
      geom_point()+
      geom_errorbarh(aes(xmax = CI95upper, xmin = CI95lower, height =0.1))+
      geom_vline(xintercept = 0, linetype='dotted', color='grey', size=1.5)+
      scale_color_manual(values=colors_list)+
      theme_classic()+
      theme(axis.text.y = element_text(size=10),
            axis.text.x = element_text(size=15),
            axis.title.x = element_text(size=15),
            axis.title.y = element_blank(),
            legend.position = "bottom",
            plot.title = element_text(size=15),
            panel.border = element_rect(fill=NA, colour="black", size=3))+
      xlim(x_min_lim, x_max_lim)+
      xlab('coefficient')
    save_plot(paste0(output_folder,plotname), coeff_p, base_height=12, base_width=11)
    
    print("Done")
    
    #### Generate prevalence plot ###
    ### Access absolute pathway abundance data 
    abs_abund <- psmelt(ps_obj)
    write.csv(abs_abund, paste0(output_folder,name_string,"_abs_abund.csv"), 
              row.names = F)
    diff_feat_abs_abund <- subset(abs_abund, OTU %in% filename$feature)
    View(diff_feat_abs_abund)
    # Filter absolute pathway abundance data by variable and level of interest
    list <- c(reference, level)
    print(list)
    diff_feat_abs_abund <- diff_feat_abs_abund %>%
      filter(!!sym(var) %in% list)
    
    # Re-order feature abs abundances according to coefficient value
    # This order is preserved in filename, which was previously ordered: 
    order <- data.frame(OTU = unique(filename$feature))
    # sort = F makes sure that the order is based on order_index df which is ordered by coefficient value
    abs_abund.df <- merge(order, diff_feat_abs_abund, by="OTU", sort = F)
    abs_abund.df$OTU <- as.factor(abs_abund.df$OTU)
    abs_abund.df <- abs_abund.df %>% mutate(OTU = fct_relevel(OTU, filename$feature))
    
    # For each significant feature in each sample, set AbunNotZero == 1 if its 
    # abundance is greater than 0. 
    abs_abund.df$AbunNotZero <- 1*(abs_abund.df$Abundance > 0)
    
    # Summarize pathway prevalence by variable of interest
    abs_abund.df.n0 <- abs_abund.df %>% 
      group_by(OTU, !!sym(var)) %>%
      summarise(N = n(), SumNotZero = sum(AbunNotZero)) %>%
      mutate(PercentNotZero = 100*(SumNotZero / N)) %>%
      ungroup()
    
    print(paste0(output_folder, "_", name_string))
    write.csv(abs_abund.df.n0, paste0(output_folder,"_",name_string,"_prevalence.csv"), 
              row.names = F)
    
    # Plot prevalence by variable of interest
    print("Plotting prevalence plot")
    
    plotname2 <- paste0(name_string,"_prev_plot.png")
    
    prev_p <- ggplot(abs_abund.df.n0, aes(x = !!sym(var),y=OTU, fill=PercentNotZero))+
      geom_tile()+
      scale_fill_distiller(palette ='RdBu', type='seq', breaks = c(0,50,100), 'Prevalence (%)')+
      theme_classic()+
      # theme(axis.text.y = element_text(size=10),
      #       axis.text.x = element_text(size=10, angle=45, hjust=1),
      #       axis.title.x = element_blank(),
      #       axis.title.y = element_blank(),
      #       plot.title = element_text(size=15),
      #       panel.border = element_rect(fill=NA, colour="black", size=3))+
      #scale_y_discrete(limits=rev)
    save_plot(paste0(output_folder,plotname2), prev_p, base_height=12, base_width=11)
    
    print("Done")
  }
}


plot_heatmap <- function(ps_obj, path_to_maaslin_hits_txt, values=c("healthy", "preclinical", "symptomatic"), color_vector=c("#444E7E","#FEB424", "#D8511D"), seed=42, output=output_path){
  ### IMPORTANT - Maaslin outputs will look like k__Bacteria.p__Firmicutes|
  ### you need to edit to it look like k__Bacteria|p__Firmicutes| - just commmand f . and replace with nothing
  set.seed(seed)
  
  
  taxa_list <- readLines(path_to_maaslin_hits_txt)
  
  temp_ps <- ps_obj %>% microViz::ps_filter(AD_Status %in% values)
  
  
  png(file=paste0(output_path,"_",deparse(substitute(ps_obj)), ".png"), width=7,height=7,units="in",res=1800)
  
  
  
  ht <-comp_heatmap(data=temp_ps, 
                    taxa=taxa_list, 
                    taxon_renamer = function(x) stringr::str_remove(x, ".*s__"),
                    sample_anno = sampleAnnotation(AD_Status = anno_sample_cat("AD_Status", legend_title = "AD_Status"))
                    
                    
  )
  draw(ht)
  dev.off()
  
}



maaslin_plots <- function(maaslin_output, ps_obj, values=c("healthy", "preclinical", "symptomatic"), metadata_val="AD_Status", n=20, num_samples=188,  color_vector=c("#444E7E", "#FEB424", "#D8511D"), output=output_path) {
  # Create output directory if it doesn't exist
  if (!file.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  sample_cutoff <- num_samples*.15
  temp_ps <- ps_obj %>% microViz::ps_filter(AD_Status %in% values)
  maaslin_df <- read_tsv(maaslin_output) %>% filter(metadata==metadata_val) %>% mutate(abs_coef=abs(coef)) %>% mutate(qval = p.adjust(pval, method = 'BH')) %>% filter(N.not.0>sample_cutoff) %>% mutate(feature=gsub(".*s__", "", feature)) %>% mutate(color_value= case_when(
    coef < 0 ~ "healthy", 
    TRUE~ value)) %>% mutate(value=as.factor(value)) %>% mutate(as.factor(color_value))
  
  write.csv(maaslin_df, paste0(output_path,"_maaslin_df_test.csv"), 
            row.names = F)
  median <- median(maaslin_df$abs_coef)
  
  
  top_hits <- maaslin_df %>%
    arrange(coef) %>%
    #mutate(pathways = fct_reorder(feature, coef)) %>%
    slice(1:n) %>% filter(coef > median)
  
  bottom_hits <- maaslin_df %>%
    arrange(coef) %>%
    #mutate(pathways = fct_reorder(feature, coef)) %>%
    slice(1:n) %>% filter(abs(coef) > median)
  
  results_df <- rbind(top_hits, bottom_hits) 
  print(dim(results_df))
  coeff_p <- ggplot(results_df, aes(x=coef, y=feature,color=color_value))+
    geom_point()+
    geom_vline(xintercept = 0, linetype='dotted', color='grey', size=.5)+
    scale_color_manual(values=color_vector)+
    theme_classic()+
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=15),
          axis.title.x = element_text(size=15),
          axis.title.y = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", size=3))+ xlab('coefficient')+facet_wrap(~ value)
  
  ggsave(filename=paste0(output_path, "_Maaslin_",deparse(substitute(ps_obj)),"_2.png"), plot=coeff_p, dpi=900, height=7, width=7)
  
  
  
}