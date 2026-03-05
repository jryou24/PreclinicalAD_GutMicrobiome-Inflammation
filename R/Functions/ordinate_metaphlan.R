ordinate.mphln <- function(phyloseq.obj, method, distance, adonisformula = NULL, 
                               group, seed=42, color_vector, markershape, saveplot,
                           output_path) {
  #ARG HOUSEKEEPING
  metadata <- data.frame(phyloseq::sample_data(phyloseq.obj))
  
  group.name <- deparse(substitute(group))
  group.idx <- grep(group.name, colnames(metadata))
  
  
  #ORDINATE
  ps.ordination <- phyloseq::ordinate(phyloseq.obj, method=method, distance=distance)
  
  #PLOT
  ordination.plot <- phyloseq::plot_ordination(phyloseq.obj, ps.ordination,
                                               type="samples", color=group.name)+
    geom_point(size=3, shape=markershape)+
    stat_ellipse()+
    scale_color_manual(values=color_vector)+ 
    theme_classic()+
    theme(axis.text.y = element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          legend.text = element_text(size=20),
          legend.position = 'left',
          plot.title = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", size=3))
  
  ordination.plot2 <- ggMarginal(ordination.plot, type='boxplot', groupFill=TRUE, 
                                 size=10)
  print("Done with PCoA")
  #SAVE PLOT
  phyloseq.obj.name <- deparse(substitute(phyloseq.obj))
  filepath <- paste(output_path, "/",
                    Sys.Date(),"_", 
                    phyloseq.obj.name,"_",
                    group.name,"_",
                    method,"_",
                    distance, ".pdf", sep="")
  
  if (saveplot == TRUE){
    ggsave(filepath, 
           plot = ordination.plot2, 
           device ='pdf', 
           width=20, height=13, units='cm')
  }
  
  #SAMPLE COORDINATES: MARGINAL t (or anova)-tests 
  coordinates <- data.frame(ps.ordination$vectors)
  print('Vector lengths: Coord axis 1, Coord axis 2, metadata$group')
  print(c(length(coordinates$Axis.1), length(coordinates$Axis.2), length(metadata[, group.idx])))
  print('Variable compared')
  print(group.name)
  
  df <- data.frame(Coord1 = coordinates$Axis.1,
                   Coord2 = coordinates$Axis.2,
                   Group = metadata[, group.idx])
  
  axis1.test <- aov(Coord1 ~ Group, data=df)
  axis2.test <- aov(Coord2 ~ Group, data=df)
  
  #ADONIS2 (PERMANOVA) TESTs
  # ps.dist <- phyloseq::distance(phyloseq.obj, method=distance) 
  # 
  # ps.adonis <- vegan::adonis2(formula(adonisformula),  
  #                             data=metadata,
  #                             na=na.omit,
  #                             permutations = 10000,
  #                             subset=complete.cases(ps.ordination))  
  
  #RETURN 
  list(ordination = ps.ordination, plot = ordination.plot2, 
       axis1test = axis1.test, axis2test = axis2.test, axistest_df = df)
  
}

ordinate.mphln.long <- function(phyloseq.obj, method, distance, adonisformula = NULL, 
                           group, seed=42, color_vector, saveplot,
                           output_path) {
  #ARG HOUSEKEEPING
  metadata <- data.frame(phyloseq::sample_data(phyloseq.obj))
  
  group.name <- deparse(substitute(group))
  group.idx <- grep(group.name, colnames(metadata))
  
  
  #ORDINATE
  ps.ordination <- phyloseq::ordinate(phyloseq.obj, method=method, distance=distance)
  
  #PLOT
  ordination.plot <- phyloseq::plot_ordination(phyloseq.obj, ps.ordination,
                                               type="samples", color=group.name)+
    geom_point(aes(shape = Time_point), size=3)+
    stat_ellipse()+
    scale_color_manual(values=color_vector)+ 
    theme_classic()+
    theme(axis.text.y = element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          legend.text = element_text(size=20),
          legend.position = 'left',
          plot.title = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", size=3))
  
  ordination.plot2 <- ggMarginal(ordination.plot, type='boxplot', groupFill=TRUE, 
                                 size=10)
  print("Done with PCoA")
  #SAVE PLOT
  phyloseq.obj.name <- deparse(substitute(phyloseq.obj))
  filepath <- paste(output_path, "/",
                    Sys.Date(),"_", 
                    phyloseq.obj.name,"_",
                    group.name,"_",
                    method,"_",
                    distance, ".pdf", sep="")
  
  if (saveplot == TRUE){
    ggsave(filepath, 
           plot = ordination.plot2, 
           device ='pdf', 
           width=20, height=13, units='cm')
  }
  
  #SAMPLE COORDINATES: MARGINAL t (or anova)-tests 
  coordinates <- data.frame(ps.ordination$vectors)
  print('Vector lengths: Coord axis 1, Coord axis 2, metadata$group')
  print(c(length(coordinates$Axis.1), length(coordinates$Axis.2), length(metadata[, group.idx])))
  print('Variable compared')
  print(group.name)
  
  df <- data.frame(Coord1 = coordinates$Axis.1,
                   Coord2 = coordinates$Axis.2,
                   Group = metadata[, group.idx])
  
  axis1.test <- aov(Coord1 ~ Group, data=df)
  axis2.test <- aov(Coord2 ~ Group, data=df)
  
  #RETURN 
  list(ordination = ps.ordination, plot = ordination.plot2, 
       axis1test = axis1.test, axis2test = axis2.test, axistest_df = df)
  
}

#' PERMANOVA for repeated measures using adonis2
#'
#' @param D An N-by-N distance matrix (must be a \code{dist} object)
#' @param subject A per-sample (length-N) character vector containing the subject identifiers
#' @param subject_data Data frame with per-subject metadata. Must have rownames matching unique(subject)
#' @param sample_data Data frame with per-sample metadata. Must have rownames matching D
#' @param metadata_order Character vector specifying the variables to include in the model
#' @param permutations Number of permutations to perform
#' @param ncores Number of cores for parallelization
#' @param by_type Type of test for adonis2 ("terms" or "margin")
#' @return An object similar to adonis2 output with empirical p-values
PERMANOVA_repeat_measures <- function(
    D,
    subject, subject_data = NULL,
    sample_data = NULL,
    metadata_order = c(names(subject_data), names(sample_data)),
    permutations = 999, ncores = 1,
    by_type = "terms") {
  
  if (!inherits(D, "dist")) {
    stop("D must be a dist object")
  }
  
  if (!is.character(subject)) stop("subject must be a character vector")
  if (length(subject) != nrow(as.matrix(D))) {
    stop("Length of subject must match number of samples in D")
  }
  if (!identical(names(subject), rownames(as.matrix(D)))) {
    stop("Names of subject vector must match sample names in D")
  }
  
  if (is.null(subject_data) & is.null(sample_data)) {
    stop("At least one of subject_data or sample_data must be provided")
  }
  
  if (is.null(subject_data)) {
    subject_data <- data.frame(.placeholder = 1, row.names = unique(subject))
  }
  
  if (length(unique(subject)) != nrow(subject_data)) {
    stop("Number of unique subjects must match rows in subject_data")
  }
  if (!setequal(unique(subject), rownames(subject_data))) {
    stop("Row names of subject_data must match unique values in subject")
  }
  
  if (is.null(sample_data)) {
    sample_data <- data.frame(.placeholder = 1, row.names = rownames(as.matrix(D)))
  }
  if (nrow(as.matrix(D)) != nrow(sample_data)) {
    stop("sample_data must have a row for each sample in D")
  }
  if (!identical(rownames(as.matrix(D)), rownames(sample_data))) {
    stop("Row names of sample_data must match row names of D")
  }
  
  if (length(intersect(names(subject_data), names(sample_data))) > 0) {
    stop("Column names must not overlap between subject_data and sample_data")
  }
  
  if (!all(metadata_order %in% c(names(subject_data), names(sample_data)))) {
    stop("metadata_order must only include columns from subject_data or sample_data")
  }
  
  mtdat <- cbind(subject_data[subject, , drop = FALSE], sample_data)
  
  if (any(is.na(mtdat[, metadata_order, drop = FALSE]))) {
    stop("Missing values found in model metadata")
  }
  
  # Construct formula explicitly
  formula <- as.formula(paste("D ~", paste(metadata_order, collapse = " + ")))
  ad <- vegan::adonis2(formula, permutations = 0, data = mtdat[, metadata_order, drop = FALSE], by = by_type)
  
  R2 <- ad$R2
  names(R2) <- rownames(ad)
  
  library(permute)
  doParallel::registerDoParallel(ncores)
  
  nullsamples <- foreach::`%dopar%`(
    foreach::foreach(i = seq_len(permutations), .combine = cbind),
    {
      subject.i <- sample(unique(subject))
      sample.i <- shuffle(nrow(sample_data), control = how(blocks = subject))
      
      i.subject_data <- subject_data[subject.i, , drop = FALSE]
      rownames(i.subject_data) <- rownames(subject_data)
      
      mtdat_perm <- cbind(i.subject_data[subject, , drop = FALSE],
                          sample_data[sample.i, , drop = FALSE])
      perm_ad <- vegan::adonis2(formula, permutations = 0, data = mtdat_perm[, metadata_order, drop = FALSE], by = by_type)
      r2_perm <- perm_ad$R2
      matrix(r2_perm, ncol = 1, dimnames = list(names(r2_perm), NULL))
    })
  
  doParallel::stopImplicitCluster()
  
  n <- length(R2)
  R2[n - 1] <- 1 - R2[n - 1]
  nullsamples[n - 1, ] <- 1 - nullsamples[n - 1, ]
  
  exceedances <- rowSums(nullsamples > R2)
  P <- (exceedances + 1) / (permutations + 1)
  P[n] <- NA
  
  ad$`Pr(>F)` <- P
  return(ad)
}

