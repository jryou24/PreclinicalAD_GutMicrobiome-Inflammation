
### Mapping the following RF Vars#######
#AD_Status':AD_Status
#LUMIPULSE_CSF.ab42.ab40.ratio:ab42_ab40_ratio
#combined.Centiloid:combined.centiloid
#Tauopathy:Tauopathy
#LUMIPULSE_CSF_pTau:LUMIPULSE_CSF_pTau
#LUMIPULSE_CSF_tTau:LUMIPULSE_CSF_tTau
#sex:sex
#age:age
#race:race
#Years.Education:EDUC
#BMI:BMI
#Hypertension:HYPERTEN
#Diabetes:DIABETES
#APOE4.status:apoe4.carriage
#####new###########
#apoe4 dose
#DEP2YRS
#ALCOHOL
#TOBAC30
#THYROID
#HYPERCHO
#CVHATT
#CBSTROKE

####Missing#############
#WMH_volume
#CortSig_Thickness
#MR_LV_RV_HIPPOCAMPUS
#Polygenic.Risk.Score
#INT_PET

#INT_MRI
#INT_CSF:Days_sinceCSF <- exclduded since im not sure about the other ones

### Using Aura's code


# FUNCTION: callBoruta helper function (called by runFeatureSelection() below)
# Required packages: Boruta, stringr
# Arguments:
#  taxa.data (dataframe)  = taxa.train.wclass (tax abundances for training cohort)
#  seed.Boruta (int)      = will be passed iteratively in defined range
#
# Return:
#  list of names of feature-selected taxa for current iteration

callBoruta <- function(taxa.data, seed.Boruta){
  set.seed(seed.Boruta)
  
  taxa.boruta <- Boruta(AD_Status~., data = taxa.data, maxRuns=500, doTrace=0)
  taxa.boruta.fix <- TentativeRoughFix(taxa.boruta)
  
  taxa.boruta.df <- data.frame('boruta' = taxa.boruta.fix$finalDecision)
  
  taxanames <- str_replace_all(rownames(subset(taxa.boruta.df, 
                                               boruta=='Confirmed')), "`", "")
  
}


# FUNCTION: runFeatureSelection - iterative feature selection function
# Required packages: Boruta, stringr
# Arguments:
#  taxa.train.data (dataframe)   = taxa.train.wclass (passed to callBoruta())
#  seed.range (num list)         = range for iteration, e.g. 1:100.
#
# Return:
#  dataframe summarizing frequency at which unique taxa were feature-selected
#  across all iterations (random seeds).

runFeatureSelection <- function(taxa.train.data, seed.range){
  
  fs.taxa.train <- vector('list', length(seed.range))
  for (i in seed.range){
    fs.taxa.train[[i]] <- callBoruta(taxa.train.data, i)
  }
  
  fs.taxa.train.summ <- data.frame(table(unlist(fs.taxa.train)))
} 


## Create data subsets (i.e. omit biomarker categories)-------------------------

# FUNCTION: createSubsets - supply clinical metadata features to omit, then 
# merge remaining clinical metadata features with previously selected taxonomic 
# features.
# Required packages: none
# Arguments:
#  RFmetadata (dataframe)    = meta.caret.im (imputed metadata)
#  nullvars (chr list)       = list of vars to exclude from the model
#  Taxa (num matrix)         = taxa.caret.boruta (abundances of selected taxa)
#
# Return:
#  Named list: Base = values for selected metadata features, WithTaxa = values
#  for selected metadata features as well as feature-selected taxa. This will
#  enable testing of improvements in model performances with addition of
#  taxonomic feature data.

createSubsets <- function(RFmetadata, nullvars, Taxa){
  if (!all(nullvars %in% colnames(RFmetadata))){
    stop("At least one variable name not in dataframe")
  }
  
  #Non-taxonomic features only  
  metadata <- RFmetadata
  metadata[ , nullvars] <- list(NULL)
  
  #With taxonomic features
  metadata.wtax <- merge(metadata, Taxa, all = TRUE, by='row.names')
  rownames(metadata.wtax) <- metadata.wtax$Row.names
  metadata.wtax$Row.names <- NULL
  
  #RETURN
  all.out <- list("Base"= metadata, 
                  "WithTaxa"= metadata.wtax)
}


train_rf_models <- function(data, control.harness, data.name.string, 
                            varsNot2Norm, shuffle.class) {
  #Cross Validation (within training cohort) 
  out <- list()
  varimportance <- list() 
  
  #Validation Set 
  out.val <- list()
  
  #Separate into train/test and validation subsets, using same index as before.
  set.seed(42)
  data1_idx <- createDataPartition(data$AD_Status, p = 0.6, list=FALSE)
  data1 <- data[data1_idx, ]
  data.val <- data[-data1_idx, ]
  
  for (i in 1:100) {
    #Create random partition of training cohort (80:20). 
    set.seed(i)
    
    train_idx <- createDataPartition(data1$AD_Status, p = 0.8, list=FALSE)
    
    data.train <- data1[train_idx, ]
    data.test <- data1[-train_idx, ]
    
    # Optionally shuffle class labels in training data (data.train).
    if (shuffle.class == TRUE) {
      data.train$AD_Status <- sample(data.train$AD_Status)
    }
    
    #Pre-process data (center and scale)
    preprocessrule <- preProcess(data.train[, !(colnames(data.train) %in% varsNot2Norm)], 
                                 method = c('center', 'scale'))
    data.train.p <- predict(preprocessrule, data.train)
    data.test.p <- predict(preprocessrule, data.test)
    
    data.val.p <- predict(preprocessrule, data.val)
    
    
    #train model
    set.seed(42)
    fit.rf <- train(AD_Status~., data = data.train.p, method = 'rf',
                    metric = 'Accuracy', trControl = control.harness)
    
    #Make predictions for this iteration's test set and store performance measures.
    predictions.rf <- predict(fit.rf, data.test.p)
    cM <- confusionMatrix(predictions.rf, data.test.p$AD_Status)
    pred.results <- c('healthy-healthy'=cM$table[1,1], 
                      'healthy-preclinical'=cM$table[2,1], 
                      'preclinical-healthy'=cM$table[1,2], 
                      'preclinical-preclinical'= cM$table[2,2], 
                      cM$overall, 
                      cM$byClass,
                      'Data' = data.name.string, 
                      'Seed' = i)
    
    out[[i]] <- pred.results
    
    #Make predictions for retained VALIDATION set and store performance measures
    predictions.val.rf <- predict(fit.rf, data.val.p)
    cM.val <- confusionMatrix(predictions.val.rf, data.val.p$AD_Status)
    pred.results.val <- c('healthy-healthy'=cM.val$table[1,1], 
                          'healthy-preclinical'=cM.val$table[2,1], 
                          'preclinical-healthy'=cM.val$table[1,2], 
                          'preclinical-preclinical'= cM.val$table[2,2], 
                          cM.val$overall, 
                          cM.val$byClass,
                          'Data' = data.name.string, 
                          'Seed' = i)
    
    out.val[[i]] <- pred.results.val
    
    
    #Find important vars
    fit.rf.importance <- varImp(fit.rf, scale=FALSE)
    fit.rf.importance.df <- fit.rf.importance$importance
    var.importance <- fit.rf.importance.df$Overall
    names(var.importance) <- row.names(fit.rf.importance.df)
    
    varimportance[[i]] <- var.importance
  }
  
  #Return 
  out.df <- data.frame(do.call('rbind', out))
  out.val.df <- data.frame(do.call('rbind', out.val))
  varimportance.df <- data.frame(do.call('rbind', varimportance))
  
  allout <- list('Pred.Results.CV'=out.df, 
                 'Pred.Results.Val'=out.val.df, 
                 'Var.Importance'=varimportance.df)
}
ml_wrapper <-  function(ps_obj, maaslin_path,
                          values=c("healthy", "preclinical"),output=output_path,
                          rf_vars=c("AD_Status", "ab42_ab40_ratio", "combined.centiloid",
                                    "Tauopathy", "LUMIPULSE_CSF_pTau", "LUMIPULSE_CSF_tTau", 
                                    "sex", "age", "race", "EDUC", "BMI", "HYPERTEN", 
                                    "DIABETES", "apoe4.carriage", "apoe.dose", "DEP2YRS",
                                    "ALCOHOL", "TOBAC30", "HYPERCHO", "THYROID", "CVHATT", "CBSTROKE", ), 
                          plot_missing_vals=TRUE, seed=42, validation_class="symptomatic"){

  
  ps_obj_name<- deparse(substitute(ps_obj))
  set.seed(42)
  temp_ps <- ps_obj %>% microViz::ps_filter(AD_Status %in% values)
  # Transform to counts, preserving precision to 0.00001%
  ps2.filt.int <- transform_sample_counts(temp_ps, function(x) trunc(x*100000))
  # Access clinical metadata
  meta.caret <- data.frame(sample_data(ps2.filt.int), check.names=FALSE)[, rf_vars]
  meta.caret[meta.caret == "NA"] <- NA
  # Access taxonomic abundance data (note transform of matrix)
  taxa.caret <- t(data.frame(otu_table(ps2.filt.int), check.names=FALSE))
  
  if(plot_missing_vals){

    meta.vimplot <- aggr(meta.caret, numbers = TRUE, sortVars = TRUE,
                       col = c("#347B98", "#66B032", "#B2D732"),
                       labels = names(meta.caret), cex.axis= 0.3, combined = FALSE)
    print(meta.vimplot)
  }
  
  meta.caret$Participant <- rownames(meta.caret)
  meta.caret.im <- kNN(meta.caret, imp_var = FALSE)
  
  
  rownames(meta.caret.im) <- meta.caret.im$Participant
  meta.caret.im$Participant <- NULL
  
  factor.vars <- c('AD_Status', 'apoe4.carriage', 'sex', 'race',
                   'HYPERTEN', 'DIABETES', "DEP2YRS",
                   "ALCOHOL", "TOBAC30", "HYPERCHO", "THYROID", "CVHATT", "CBSTROKE")
  
  meta.caret.im[, factor.vars] <- lapply(meta.caret.im[, factor.vars], factor)
  train_idx <- createDataPartition(meta.caret.im$AD_Status, p = 0.6, list=FALSE)
  
  # Access metadata (biomarkers) and taxanomic abundances for training cohort.
  meta.train <- meta.caret.im[train_idx, ]
  taxa.train <- taxa.caret[train_idx, ] 
  
  # Merge sample class identity with taxonomic abundance data 
  class <- subset(meta.train, select=c('AD_Status'))
  
  taxa.train.wclass <- merge(class, taxa.train, by='row.names')
  
  rownames(taxa.train.wclass) <- taxa.train.wclass$Row.names
  taxa.train.wclass$Row.names <- NULL
  
  
  # Carry out iterative feature selection (may take a little while..~20 min)
  fs.taxa <- runFeatureSelection(taxa.train.wclass, 1:100)
  # Recommended to save workspace at this point.
  # Filter for taxa selected in > 25% of iterations.
  fs.taxa.top <- subset(fs.taxa, Freq > 25, select='Var1')
  # Subset taxa abundance data to these selected taxa (for entire cohort including
  # train and validation sets -- will be re-partitioned later using same index). 
  taxa.caret.boruta <- taxa.caret[ , colnames(taxa.caret) %in% fs.taxa.top$Var1]
  
  # Subset taxa abundance data to these selected taxa (for entire cohort including
  # train and validation sets -- will be re-partitioned later using same index). 
  taxa.caret.boruta <- taxa.caret[ , colnames(taxa.caret) %in% fs.taxa.top$Var1]
  
  ## Plot relative abundances of feature-selected taxa----------------------------
  
  # Subset melted phyloseq dataframe to the feature-selected taxa.
  ps2.filt.int.df <- psmelt(ps2.filt.int) 
  
  diffspecies.all.df <- read_tsv(maaslin_path)
  diffspecies.a.df <- subset(diffspecies.all.df, metadata == 'AD_Status' & 
                               N.not.0 > 25 & abs(coef) > 0.15 & value=='preclinical') 
  
  diffspecies.a.df <- diffspecies.a.df %>% arrange(desc(coef))
  diffspecies.a <- unique(diffspecies.a.df$feature)

  ps2.filt.int.a.df <- subset(ps2.filt.int.df, OTU %in% diffspecies.a.df$feature)
  ps2.filt.int.boruta <- subset(ps2.filt.int.df, OTU %in% fs.taxa.top$Var1)
  
  
  # Make facet labels vector 
  # List of unique OTUs, split for just species name, but name the species vector 
  # with full OTU name.
  taxafeature.labels.df <- fs.taxa.top %>% separate(Var1, c(NA, 'species'), '\\|s__')
  taxafeature.labels <- taxafeature.labels.df$species
  names(taxafeature.labels) <- fs.taxa.top$Var1
  
  
  # Plot relative abundances of feature-selected taxa.
  # Note transformation back to percent scale, with addition of small pseudo-count
  # to enable plotting on a log scale. 
  p_taxafeature <- ggplot(ps2.filt.int.boruta, 
                          aes(x=AD_Status, y=Abundance/100000+0.0001))+
    geom_violin(aes(color=AD_Status), draw_quantiles = c(0.25,0.75))+
    geom_jitter(aes(color=AD_Status),width=0.2, size=1, alpha=0.8, pch=21)+
    scale_color_brewer(palette ='Dark2')+
    theme_classic()+
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10, angle=45, hjust=1),
          axis.title.x = element_blank(10),
          axis.title.y = element_text(size=10),
          legend.position = "none",
          plot.title = element_text(size=10),
          panel.border = element_rect(fill=NA, colour="black", size=1))+
    facet_wrap(~OTU, labeller = labeller(OTU = taxafeature.labels), nrow=1)+
    scale_y_continuous(trans="log10", labels=comma) +
    ylab('Abundance (%)')
  
  ggsave(filename=paste0(output,"_",ps_obj_name, "_taxafeature.png"), plot=p_taxafeature, dpi=900, height=7, width=14)

  
  
  }

