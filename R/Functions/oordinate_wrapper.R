pca_iris_wrapper <- function(ps_obj, 
        values=c("healthy", "preclinical", "symptomatic"), 
        color_vector=c("#444E7E","#FEB424", "#D8511D"), output=output_path){
  ps_obj_name<- deparse(substitute(ps_obj))
  temp_ps <- ps_obj %>% microViz::ps_filter(AD_Status %in% values)
  ### PCA plot
  pca <- temp_ps %>% microViz::tax_transform("identity", rank = "Species", zero_replace = 0) %>% ord_calc() %>%  ord_plot(color = "AD_Status", shape = "Time_point", plot_taxa = 1:5, size = 2) +
    scale_color_manual(values=color_vector)+theme_bw() +
    ggside::geom_xsidedensity(aes(fill = AD_Status,), alpha = 0.5, show.legend = FALSE) +
    ggside::geom_ysidedensity(aes(fill = AD_Status), alpha = 0.5, show.legend = FALSE) +
    ggside::theme_ggside_void()+stat_ellipse(aes(color = AD_Status)) +scale_fill_manual(values=color_vector)
  
  print("PCA Comparisons for ")
  print(deparse(substitute(ps_obj)))
  print(TukeyHSD(aov(PC1 ~ AD_Status, data=pca$data)))
  print(TukeyHSD(aov(PC2 ~ AD_Status, data=pca$data)))
  ### USELESS IRIS 
  print("Generating Iris")
  iris <- temp_ps %>% microViz::tax_transform("identity", rank = "Species", zero_replace =  0) %>% ord_calc() %>% ord_plot_iris(tax_level = "Species",n_taxa=10, anno_binary = "AD_Status", anno_colour = "AD_Status", anno_binary_style = list())+  theme(legend.text = element_text(size = 5))+
    scale_colour_manual(values = color_vector)
  
  pca_path <- paste(c(output, ps_obj_name, "pca.png"), collapse="_")
  iris_path <- paste(c(output, ps_obj_name, "iris.png"), collapse="_")
  
  ggsave(filename=pca_path, plot=pca, dpi=900, height=7, width=7)
  ggsave(filename=iris_path, plot=iris, dpi=900, height=7, width=7)
}

unifrac_wrapper <- function(ps_obj,
                            values=c("healthy", "preclinical", "symptomatic"), 
                            color_vector=c("#444E7E","#FEB424", "#D8511D"),
                            seed=42, output=output_path){
  
  temp_ps <- ps_obj %>% microViz::ps_filter(AD_Status %in% values)
  #### UNWEIGHTED UNIFRAC
  print("Generating Unifrac")
  pcoa_unweighted_unifrac <- temp_ps %>%tax_transform("identity", rank = "unique",  zero_replace = 0) %>%
    dist_calc("unifrac") %>% ord_calc(method = "PCoA") %>%
    ord_plot(color = "AD_Status", shape = "Time_point", size = 2) +
    scale_color_manual(values=color_vector)+theme_bw() +
    ggside::geom_xsidedensity(aes(fill = AD_Status,), alpha = 0.5, show.legend = FALSE) +
    ggside::geom_ysidedensity(aes(fill = AD_Status), alpha = 0.5, show.legend = FALSE) +
    ggside::theme_ggside_void()+stat_ellipse(aes(color = AD_Status)) + scale_fill_manual(values=color_vector)
  print("Unweighted Unifrac PCoA Comparisons for ")
  print(deparse(substitute(ps_obj)))
  print(TukeyHSD(aov(MDS1 ~ AD_Status, data=pcoa_unweighted_unifrac$data)))
  print(TukeyHSD(aov(MDS2 ~ AD_Status, data=pcoa_unweighted_unifrac$data)))
  
  pcoa_path <- paste(c(output, deparse(substitute(ps_obj)), "UNIFRAC_pcoa.png"), collapse="_")
  ggsave(filename=pcoa_path, plot=pcoa_unweighted_unifrac, dpi=900, height=7, width=7)
}



bray_curtis_wrapper <- function(ps_obj,
                            values=c("healthy", "preclinical", "symptomatic"), 
                            color_vector=c("#444E7E","#FEB424", "#D8511D"),
                            seed=42, output=output_path){
  
  temp_ps <- ps_obj %>% microViz::ps_filter(AD_Status %in% values)
  
  pcoa_bc <- temp_ps %>%tax_transform("identity", rank = "unique",  zero_replace = 0) %>%
    dist_calc("bray") %>% ord_calc(method = "PCoA") %>%
  ord_plot(color = "AD_Status", shape = "Time_point", size = 2) +
    scale_color_manual(values=color_vector)+theme_bw() +
    ggside::geom_xsidedensity(aes(fill = AD_Status,), alpha = 0.5, show.legend = FALSE) +
    ggside::geom_ysidedensity(aes(fill = AD_Status), alpha = 0.5, show.legend = FALSE) +
    ggside::theme_ggside_void()+stat_ellipse(aes(color = AD_Status)) + scale_fill_manual(values=color_vector)
  print("Bray-Curtis PCoA Comparisons for ")
  print(deparse(substitute(ps_obj)))
  print(TukeyHSD(aov(MDS1 ~ AD_Status, data=pcoa_bc$data)))
  print(TukeyHSD(aov(MDS2 ~ AD_Status, data=pcoa_bc$data)))

  pcoa_path <- paste(c(output, deparse(substitute(ps_obj)), "BC_pcoa.png"), collapse="_")
  ggsave(filename=pcoa_path, plot=pcoa_bc, dpi=900, height=7, width=7)
}

alpha_div_wrapper <- function(ps_obj,
                              values=c("healthy", "preclinical", "symptomatic"),
                              save_adiv_plot = FALSE,
                              color_vector=c("healthy" = "#49b9a5", "preclinical" = "#eebe61", "symptomatic" = "#a98ec2"),
                              seed=42, output=output_path)
{
  temp_ps <- ps_obj %>% microViz::ps_filter(AD_Status %in% values)
  #### Alpha Diversity plots
  ps_obj_transformed <- transform_sample_counts(temp_ps, function(x) trunc(x*100000))
  adiv <- estimate_richness(ps_obj_transformed, measures=c('Observed', 'Shannon'))
  
  sample_data(temp_ps)$Richness <- adiv$Observed
  sample_data(temp_ps)$Shannon <- adiv$Shannon

  if(save_adiv_plot == TRUE){
    ps2_obj_df <- data.frame(sample_data(temp_ps))
    ps2_obj_df_gathered <- gather(ps2_obj_df, Alpha_Measure, Value, Richness:Shannon, factor_key=TRUE)
    
    alpha_div <- ggplot(ps2_obj_df_gathered, aes(x=AD_Status, y=Value, color=AD_Status))+
      geom_boxplot(outlier.shape=NA)+
      geom_jitter(width=0.2)+
      scale_color_manual(values=color_vector)+
      facet_wrap(~Alpha_Measure, nrow=1, scales = 'free_y')+
      theme_classic()+
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=20),
            axis.text.y = element_text(size=20),
            axis.text.x = element_text(size=15, angle=45, hjust=1),
            strip.text.x = element_text(size=20),
            panel.border = element_rect(fill=NA, colour="black", size=2))+
      ylab('Taxa (MetaPhlAn3)')+
      scale_y_continuous(limits=c(0, NA))
    print("Alpha Div")
    print(deparse(substitute(ps_obj)))
    print(TukeyHSD(aov(Richness ~ AD_Status, data=ps2_obj_df)))
    print(TukeyHSD(aov(Shannon ~ AD_Status, data=ps2_obj_df)))
    
    alpha_path <- paste(c(output, deparse(substitute(ps_obj)), "alpha_div.png"), collapse="_")
    ggsave(filename=alpha_path, plot=alpha_div, dpi=900, height=7, width=9)
    
  }
 
  return(temp_ps)
  
}


alpha_div_wrapper_adjust_var <- function(ps_obj,
                                         variable="AD_Status",
                                         save_adiv_plot = FALSE,
                                         values=c("healthy", "preclinical", "symptomatic"), 
                                         color_vector=c("healthy" = "#49b9a5", "preclinical" = "#eebe61", "symptomatic" = "#a98ec2"),
                                         seed=42, output=output_path) {
  
  temp_ps <- ps_obj %>% microViz::ps_filter(!!sym(variable) %in% values)
  #### Alpha Diversity plots
  ps_obj_transformed <- transform_sample_counts(temp_ps, function(x) trunc(x*100000))
  adiv <- estimate_richness(ps_obj_transformed, measures=c('Observed', 'Shannon','Simpson'))
  
  sample_data(temp_ps)$Richness <- adiv$Observed
  sample_data(temp_ps)$Shannon <- adiv$Shannon
  sample_data(temp_ps)$Simpson <- adiv$Simpson
  
  if(save_adiv_plot == TRUE){
    ps2_obj_df <- data.frame(sample_data(temp_ps))
    ps2_obj_df_gathered <- gather(ps2_obj_df, Alpha_Measure, Value, Richness:Simpson, factor_key=TRUE)
    
    alpha_div <- ggplot(ps2_obj_df_gathered, aes(x=!!sym(variable), y=Value, color=!!sym(variable)))+
      geom_boxplot(outlier.shape=NA)+
      geom_jitter(width=0.2)+
      scale_color_manual(values=color_vector)+
      facet_wrap(~Alpha_Measure, nrow=1, scales = 'free_y')+
      theme_classic()+
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=20),
            axis.text.y = element_text(size=20),
            axis.text.x = element_text(size=15, angle=45, hjust=1),
            strip.text.x = element_text(size=20),
            panel.border = element_rect(fill=NA, colour="black", size=2))+
      ylab('Taxa (MetaPhlAn3)')+
      scale_y_continuous(limits=c(0, NA))
    alpha_path <- paste(c(output, deparse(substitute(ps_obj)),variable,"alpha_div.pdf"), collapse="_")
    ggsave(filename=alpha_path, plot=alpha_div, bg= "transparent", height=7, width=7)

    print("Alpha Div")
    print(deparse(substitute(ps_obj)))
    # Dynamically create the formula for Richness
    richness_formula <- as.formula(paste("Richness ~", variable))
    # Perform the ANOVA and Tukey's HSD test for Richness
    richness_aov <- aov(richness_formula, data = ps2_obj_df)
    print(TukeyHSD(richness_aov))
    # Dynamically create the formula for Shannon
    shannon_formula <- as.formula(paste("Shannon ~", variable))
    # Perform the ANOVA and Tukey's HSD test for Shannon
    shannon_aov <- aov(shannon_formula, data = ps2_obj_df)
    print(TukeyHSD(shannon_aov))
    # Dynamically create the formula for Simpson
    simpson_formula <- as.formula(paste("Simpson ~", variable))
    # Perform the ANOVA and Tukey's HSD test for Simpson
    simpson_aov <- aov(simpson_formula, data = ps2_obj_df)
    print(TukeyHSD(simpson_aov))
    
  }
  
  return(ps2_obj_df)
  
}

unifrac_wrapper_adjust_var <- function(ps_obj,
                            values=c("healthy", "preclinical", "symptomatic"), 
                            color_vector=c("#444E7E","#FEB424", "#D8511D"),
                            var_for_color = "AD_Status",
                            seed=42, output=output_path){
  
  temp_ps <- ps_obj %>% microViz::ps_filter(!!sym(var_for_color) %in% values)
  #### UNWEIGHTED UNIFRAC
  print("Generating Unifrac")
  pcoa_unweighted_unifrac <- temp_ps %>%tax_transform("identity", rank = "unique",  zero_replace = 0) %>%
    dist_calc("unifrac") %>% ord_calc(method = "PCoA") %>%
    ord_plot(color = var_for_color, size = 2) +
    scale_color_manual(values=color_vector)+theme_bw() +
    ggside::geom_xsidedensity(aes(fill = !!sym(var_for_color),), alpha = 0.5, show.legend = FALSE) +
    ggside::geom_ysidedensity(aes(fill = !!sym(var_for_color)), alpha = 0.5, show.legend = FALSE) +
    ggside::theme_ggside_void()+stat_ellipse(aes(color = !!sym(var_for_color))) + scale_fill_manual(values=color_vector)
  print("Unweighted Unifrac PCoA Comparisons for ")
  print(deparse(substitute(ps_obj)))
  print(TukeyHSD(aov(MDS1 ~ AD_Status, data=pcoa_unweighted_unifrac$data)))
  print(TukeyHSD(aov(MDS2 ~ AD_Status, data=pcoa_unweighted_unifrac$data)))
  
  pcoa_path <- paste(c(output, deparse(substitute(ps_obj)),var_for_color,"UNIFRAC_pcoa.png"), collapse="_")
  ggsave(filename=pcoa_path, plot=pcoa_unweighted_unifrac, dpi=900, height=7, width=7)
  
  return(pcoa_unweighted_unifrac)
}




permanova_wrapper <- function(ps_obj, group,values=c("healthy", "preclinical", "symptomatic"), var = "AD_Status",
            color_vector=c("#444E7E","#FEB424", "#D8511D"), dist="unifrac", ss_method = "margin",
            seed=42, output=output_path, n_perm=10000, permnova=TRUE, 
            impute=0, adonis_values=c("sex", "age", "apoe", "race", "EDUC", "DIABETES", "BMI", "cancer", "stool_date_cont", "adi_natrank2019", "mrfei", "AD_Status"))
{
  temp_ps <- ps_obj %>% microViz::ps_filter(!!sym(var) %in% values)
  permanova<- temp_ps %>%tax_transform("identity", rank = "unique",  zero_replace =  impute) %>%
    dist_calc(dist) %>%
    ord_calc(method = "PCoA") %>%
    dist_permanova(
      seed = seed, # for set.seed to ensure reproducibility of random process
      n_processes = 8, n_perms = n_perm, # you should use at least 999!
      variables = adonis_values,
      by = ss_method
    )
  print(permanova)
}


compare_beta_div <- function(ps_obj, distance_metric, distance_metric_string, output_name, plot = T){
  # Create output directory if it doesn't exist
  output_path <- paste0(output_path,"_beta_comparisons/")
  if (!file.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  # obtain metadata table from phyloseq object
  meta <- data.frame(phyloseq::sample_data(ps_obj))
  #View(meta)
  # getting rid of all date columns because format got messed up 
  meta$Date_biomarker_measured <- NULL
  meta$Stool_Date_Collected <- NULL
  meta$BIRTH <- NULL
  meta$TESTDATE.cdr <- NULL
  meta$PET_Date.pib <- NULL
  meta$PET_Date.av45 <- NULL
  meta$Plasma_Date <- NULL
  meta$PET_Date <- NULL
  meta$TESTDATE.health <- NULL
  
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
  # Merge final bc distance matrix with metadata
  ps.dist.pairwise_uniq <- merge(ps.dist.pairwise_uniq, meta, by.x = "Sample_ID1", by.y = "Sample_ID")
  ps.dist.pairwise_uniq <- merge(ps.dist.pairwise_uniq, meta, by.x = "Sample_ID2", by.y = "Sample_ID")
  
  ### Comparing within and between participant distances across AD_Status groups 
  # Adding comparison type column 
  ps.dist.pairwise_uniq$Comparison[ps.dist.pairwise_uniq$ID.x == ps.dist.pairwise_uniq$ID.y] <- "Within participant"
  ps.dist.pairwise_uniq$Comparison[ps.dist.pairwise_uniq$ID.x != ps.dist.pairwise_uniq$ID.y] <- "Between participant"
  ps.dist.pairwise_uniq$AD_Status_group[ps.dist.pairwise_uniq$AD_Status.x == ps.dist.pairwise_uniq$AD_Status.y & ps.dist.pairwise_uniq$AD_Status.x == "healthy"] <- "healthy"
  ps.dist.pairwise_uniq$AD_Status_group[ps.dist.pairwise_uniq$AD_Status.x == ps.dist.pairwise_uniq$AD_Status.y & ps.dist.pairwise_uniq$AD_Status.x == "preclinical"] <- "preclinical"
  ps.dist.pairwise_uniq$AD_Status_group[ps.dist.pairwise_uniq$AD_Status.x == ps.dist.pairwise_uniq$AD_Status.y & ps.dist.pairwise_uniq$AD_Status.x == "symptomatic"] <- "symptomatic"
  ps.dist.pairwise_uniq$AD_Status_group[ps.dist.pairwise_uniq$AD_Status.x == ps.dist.pairwise_uniq$AD_Status.y & ps.dist.pairwise_uniq$AD_Status.x == "other"] <- "other"
  ps.dist.pairwise_uniq$AD_Status_group[ps.dist.pairwise_uniq$AD_Status.x != ps.dist.pairwise_uniq$AD_Status.y] <- "different"
  ps.dist.pairwise_uniq$AD_Status_group <- factor(ps.dist.pairwise_uniq$AD_Status_group, levels = c("healthy", "preclinical", "symptomatic","other","different"))
  write.csv(ps.dist.pairwise_uniq, paste0(output_path, deparse(substitute(ps_obj)), "_", distance_metric_string, "_",
                                          output_name, ".csv"))
  #View(ps.dist.pairwise_uniq)
  #### PLOT #####
  if(plot==TRUE){
    
    comparison_plot <- ggplot(data = ps.dist.pairwise_uniq, aes(x = Comparison, y = !!sym(distance_metric), color = AD_Status_group)) + 
      geom_boxplot() + 
      scale_x_discrete(guide = guide_axis(angle = 60)) +
      scale_color_manual(values=natparks.pals("Yellowstone")) +
      facet_wrap(~ps.dist.pairwise_uniq$AD_Status_group, nrow = 1) + 
      ylab(distance_metric_string) + theme_test() +
      ylim(0,1.2) +
      labs(title = paste0(distance_metric_string, " within and between participants across AD Status Groups"))
    ggsave(paste0(output_path,output_name,"_",distance_metric,"_compare_AD_Status.pdf"), plot = comparison_plot, bg = "transparent", width = 10, height = 8) 
  }
  
  ## Prepare for Wilcoxon ranksum tests between AD Status and participant pair groups
  ps.dist.pairwise_uniq$ADgroup_pair <- paste(ps.dist.pairwise_uniq$AD_Status_group, ps.dist.pairwise_uniq$Comparison)
  ps.dist.pairwise_uniq$ADgroup_pair <- factor(ps.dist.pairwise_uniq$ADgroup_pair, levels = c("healthy Between participant", "healthy Within participant", "preclinical Between participant", "preclinical Within participant",
                                                                                              "symptomatic Between participant", "symptomatic Within participant",
                                                                                              "other Between participant","other Within participant",
                                                                                              "different Between participant", "different Within participant"))
  
  return(ps.dist.pairwise_uniq)
}




##Based on ordinate2.mphlan.AF 
ordinate.wrapper <- function(ps_obj, group,values=c("healthy", "preclinical", "symptomatic"), color_vector=c("#444E7E","#FEB424", "#D8511D"),seed=42, output=output_path, n_perm=8, permnova=TRUE, impute=0, adonis_values=c("sex", "age", "apoe", "race", "EDUC", "DIABETES", "BMI", "cancer", "stool_date_cont", "ADI_STATERNK_19", "mrfei", "AD_Status"))
{
  ps_obj_name<- deparse(substitute(ps_obj))
  # Filter our to only include AD status we care about 
  ## Skipping tax_filter since that seems to be focused on raw read counts
  #rank = unique which does not aggregate taxa, since we already filtered out other taxa above and below species level we should be good 
  temp_ps <- ps_obj %>% microViz::ps_filter(AD_Status %in% values)
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ### PCA plot
  pca <- temp_ps %>% microViz::tax_transform("identity", rank = "Species", zero_replace = impute) %>% ord_calc() %>%  ord_plot(color = "AD_Status", shape = "Time_point", plot_taxa = 1:5, size = 2) +
    scale_color_manual(values=color_vector)+theme_bw() +
    ggside::geom_xsidedensity(aes(fill = AD_Status,), alpha = 0.5, show.legend = FALSE) +
    ggside::geom_ysidedensity(aes(fill = AD_Status), alpha = 0.5, show.legend = FALSE) +
    ggside::theme_ggside_void()+stat_ellipse(aes(color = AD_Status)) +scale_fill_manual(values=color_vector)
  
  print("PCA Comparisons for ")
  print(deparse(substitute(ps_obj)))
  print(TukeyHSD(aov(PC1 ~ AD_Status, data=pca$data)))
  print(TukeyHSD(aov(PC2 ~ AD_Status, data=pca$data)))
  ### USELESS IRIS 
  print("Generating Iris")
  iris <- temp_ps %>% microViz::tax_transform("identity", rank = "Species", zero_replace =  impute) %>% ord_calc() %>% ord_plot_iris(tax_level = "Species",n_taxa=10, anno_binary = "AD_Status", anno_colour = "AD_Status", anno_binary_style = list())+  theme(legend.text = element_text(size = 5))+
    scale_colour_manual(values = color_vector)
  
  #      print("Generating Barplot")
  bar_plot <- temp_ps  %>% microViz::tax_transform("identity", rank = "Species", zero_replace =  impute) %>%
    comp_barplot("Species", n_taxa = 15, merge_other = FALSE, label = NULL, tax_order=prev) +
    facet_wrap(vars(AD_Status), scales = "free") + # scales = "free" is IMPORTANT!
    coord_flip() + ggtitle(
      "Species composition by AD Status " ) +theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))  
  
  
  
  
  
  #### UNWEIGHTED UNIFRAC
  print("Generating Unifrac")
  pcoa_unweighted_unifrac <- temp_ps %>%tax_transform("identity", rank = "unique",  zero_replace =  impute) %>%
    dist_calc("unifrac") %>%
    ord_calc(method = "PCoA") %>%
    ord_plot(color = "AD_Status", shape = "Time_point", plot_taxa = 1:5, size = 2) +
    scale_color_manual(values=color_vector)+theme_bw() +
    ggside::geom_xsidedensity(aes(fill = AD_Status,), alpha = 0.5, show.legend = FALSE) +
    ggside::geom_ysidedensity(aes(fill = AD_Status), alpha = 0.5, show.legend = FALSE) +
    ggside::theme_ggside_void()+stat_ellipse(aes(color = AD_Status)) + scale_fill_manual(values=color_vector)
  
  pcoa_unweighted_unifrac_barplot <- temp_ps %>%tax_transform("identity", rank = "unique",  zero_replace =  impute) %>%
    dist_calc("unifrac") %>%
    ord_calc(method = "PCoA") %>%
    ord_plot(color = "AD_Status", shape = "Time_point", plot_taxa = 1:5, size = 2) +
    scale_color_manual(values=color_vector)+theme_bw() +scale_fill_manual(values=color_vector)+
    ggside::geom_xsideboxplot(aes(fill = AD_Status, y = AD_Status), orientation = "y") +
    ggside::geom_ysideboxplot(aes(fill = AD_Status, x = AD_Status), orientation = "x") +
    ggside::scale_xsidey_discrete(labels = NULL) +
    ggside::scale_ysidex_discrete(labels = NULL) +
    ggside::theme_ggside_void()+stat_ellipse(aes(color = AD_Status))
  
  print("Unweighted Unifrac PCoA Comparisons for ")
  print(deparse(substitute(ps_obj)))
  print(TukeyHSD(aov(MDS1 ~ AD_Status, data=pcoa_unweighted_unifrac$data)))
  print(TukeyHSD(aov(MDS2 ~ AD_Status, data=pcoa_unweighted_unifrac$data)))
  
  
  permanova<- temp_ps %>%tax_transform("identity", rank = "unique",  zero_replace =  impute) %>%
    dist_calc("unifrac") %>%
    ord_calc(method = "PCoA") %>%
    dist_permanova(
      seed = 42, # for set.seed to ensure reproducibility of random process
      n_processes = 8, n_perms = 10000, # you should use at least 999!
      variables = adonis_values
    )
  print(permanova)
  
  #### Alpha Diversity plots
  ps_obj_transformed <- transform_sample_counts(temp_ps, function(x) trunc(x*100000))
  adiv <- estimate_richness(ps_obj_transformed, measures=c('Observed', 'Shannon'))
  
  sample_data(temp_ps)$Richness <- adiv$Observed
  sample_data(temp_ps)$Shannon <- adiv$Shannon
  
  ps2_obj_df <- data.frame(sample_data(temp_ps))
  ps2_obj_df_gathered <- gather(ps2_obj_df, Alpha_Measure, Value, Richness:Shannon, factor_key=TRUE)
  
  alpha_div <- ggplot(ps2_obj_df_gathered, aes(x=AD_Status, y=Value, color=AD_Status))+
    geom_boxplot(outlier.shape=NA)+
    geom_jitter(width=0.2)+
    scale_color_manual(values=color_vector)+theme_bw() +
    facet_wrap(~Alpha_Measure, nrow=1, scales = 'free_y')+
    theme_classic()+
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.text.x = element_text(size=15, angle=45, hjust=1),
          strip.text.x = element_text(size=20),
          panel.border = element_rect(fill=NA, colour="black", size=2))+
    ylab('Taxa (MetaPhlAn4)')+
    scale_y_continuous(limits=c(0, NA))
  print("Alpha Div")
  print(deparse(substitute(ps_obj)))
  print(TukeyHSD(aov(Richness ~ AD_Status, data=ps2_obj_df)))
  print(TukeyHSD(aov(Shannon ~ AD_Status, data=ps2_obj_df)))
  
  
  
  
  
  plot_vector <- list(pca, iris, pcoa_unweighted_unifrac, bar_plot, alpha_div, pcoa_unweighted_unifrac_barplot)
  plot_vector_names <- c("pca", "iris", "pcoa_unweighted_unifrac", "bar_plot", "alpha_div", "pcoa_unweighted_unifrac_barplot")
  
  for(i in 1:length(plot_vector)){
    print(i)
    path <- paste(c(output, ps_obj_name, plot_vector_names[i], ".png"), collapse="_")
    print(path)
    
    ggsave(filename=path, plot=plot_vector[[i]], dpi=900, height=7, width=7)
  }
}







humann.ordinate.wrapper <- function(ps_obj, adonis_values=c("sex", "age", "apoe", "race", "EDUC", "DIABETES", "BMI", "cancer", "stool_date_cont", "ADI_STATERNK_19", "mrfei", "AD_Status"), values=c("healthy", "preclinical", "symptomatic"), group, color_vector=c("#444E7E","#FEB424", "#D8511D"),seed=42, output=output_path, n_perm=8, permnova=TRUE, impute=0){
  
  ps_obj_name<- deparse(substitute(ps_obj))
  temp_ps <- ps_obj %>% microViz::ps_filter(AD_Status %in% c("healthy", "preclinical", "symptomatic"))
  print("Generating Bray Curtis")
  
  pcoa_bray <- temp_ps %>%tax_transform("identity", rank = "unique",  zero_replace =  impute) %>%
    dist_calc("bray") %>%
    ord_calc(method = "PCoA") %>%
    ord_plot(color = "AD_Status", shape = "Time_point", plot_taxa = 1:5, size = 2) +
    scale_color_manual(values=color_vector)+theme_bw() +
    ggside::geom_xsidedensity(aes(fill = AD_Status,), alpha = 0.5, show.legend = FALSE) +
    ggside::geom_ysidedensity(aes(fill = AD_Status), alpha = 0.5, show.legend = FALSE) +
    ggside::theme_ggside_void()+stat_ellipse(aes(color = AD_Status)) + scale_fill_manual(values=color_vector)
  
  print("Bray Curtis PCoA Comparisons for ")
  print(deparse(substitute(ps_obj)))
  print(TukeyHSD(aov(MDS1 ~ AD_Status, data=pcoa_bray$data)))
  print(TukeyHSD(aov(MDS2 ~ AD_Status, data=pcoa_bray$data)))
  
  
  pcoa_bray_barplot <- temp_ps %>%tax_transform("identity", rank = "unique",  zero_replace =  impute) %>%
    dist_calc("bray") %>%
    ord_calc(method = "PCoA") %>%
    ord_plot(color = "AD_Status", shape = "Time_point", plot_taxa = 1:5, size = 2) +
    scale_color_manual(values=color_vector)+theme_bw() +scale_fill_manual(values=color_vector)+
    ggside::geom_xsideboxplot(aes(fill = AD_Status, y = AD_Status), orientation = "y") +
    ggside::geom_ysideboxplot(aes(fill = AD_Status, x = AD_Status), orientation = "x") +
    ggside::scale_xsidey_discrete(labels = NULL) +
    ggside::scale_ysidex_discrete(labels = NULL) +
    ggside::theme_ggside_void()+stat_ellipse(aes(color = AD_Status))
  
  print("Bray Curtis PCoA Comparisons for ")
  print(deparse(substitute(ps_obj)))
  print(TukeyHSD(aov(MDS1 ~ AD_Status, data=pcoa_bray$data)))
  print(TukeyHSD(aov(MDS2 ~ AD_Status, data=pcoa_bray$data)))
  #### Alpha Diversity plots
  ps_obj_transformed <- transform_sample_counts(temp_ps, function(x) trunc(x*100000))
  adiv <- estimate_richness(ps_obj_transformed, measures=c('Observed', 'Shannon'))
  
  sample_data(temp_ps)$Richness <- adiv$Observed
  sample_data(temp_ps)$Shannon <- adiv$Shannon
  
  ps2_obj_df <- data.frame(sample_data(temp_ps))
  ps2_obj_df_gathered <- gather(ps2_obj_df, Alpha_Measure, Value, Richness:Shannon, factor_key=TRUE)
  
  
  
  alpha_div <- ggplot(ps2_obj_df_gathered, aes(x=AD_Status, y=Value, color=AD_Status))+
    geom_boxplot(outlier.shape=NA)+
    geom_jitter(width=0.2)+
    scale_color_manual(values=color_vector)+theme_bw() +
    facet_wrap(~Alpha_Measure, nrow=1, scales = 'free_y')+
    theme_classic()+
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.text.x = element_text(size=15, angle=45, hjust=1),
          strip.text.x = element_text(size=20),
          panel.border = element_rect(fill=NA, colour="black", size=2))+
    ylab('Pathways (Humann3)')+
    scale_y_continuous(limits=c(0, NA))
  print("Alpha Div")
  print(deparse(substitute(ps_obj)))
  print(TukeyHSD(aov(Richness ~ AD_Status, data=ps2_obj_df)))
  print(TukeyHSD(aov(Shannon ~ AD_Status, data=ps2_obj_df)))
  
  
  
  
  plot_vector <- list(pcoa_bray, alpha_div, pcoa_bray_barplot)
  plot_vector_names <- c("bray_humann",  "alpha_div_humann", "pcoa_bray_barplot")
  
  for(i in 1:length(plot_vector)){
    print(i)
    path <- paste(c(output, ps_obj_name, plot_vector_names[i], ".png"), collapse="_")
    print(path)
    
    ggsave(filename=path, plot=plot_vector[[i]], dpi=900, height=7, width=7)
  }
}

