plot_regression <- function(df, x_var, y_var, color_var, x_label, y_label, custom_colors, model0_formula,
                            model_formula, save_plot = TRUE, save_ANOVA = TRUE) {
  
  # Create plot dynamically
  plot <- ggplot(df, aes_string(x = x_var, y = y_var, color = color_var)) +
    geom_smooth(method = 'lm', show.legend = TRUE) +
    geom_point(show.legend = FALSE) +
    scale_color_manual(values = custom_colors) +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_text(size = 15),
      strip.text.x = element_text(size = 20),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 2)
    ) +
    ylim(0, NA) +
    ylab(y_label) +
    xlab(x_label)
  
  # Save plot if save_plot is TRUE
  if (save_plot) {
    filename <- paste0(Sys.Date(), "_", color_var, "_", x_var, "_", y_var, "_regression.pdf")
    ggsave(filename, plot = plot, bg = "transparent", width = 8, height = 6)
  }
  
  # Fit models
  ## Remove NAs
  df_complete <- na.omit(df[, c(x_var, y_var, color_var)])
  model0 <- lm(as.formula(model0_formula), data = df_complete, na.action = na.exclude)
  model <- lm(as.formula(model_formula), data = df_complete, na.action = na.exclude)
  # Perform ANOVA
  aov_result <- anova(model0, model)
  aov_df <- as.data.frame(aov_result)
  
  if (save_ANOVA) {
    # Save ANOVA results
    anova_filename <- paste0(Sys.Date(), "_", color_var, "_", x_var, "_", y_var, "_regression_anova.csv")
    write.csv(aov_df, anova_filename, row.names = FALSE)
  }
  return(list(plot = plot, anova = aov_df))
}

# Example usage:
# plot_regression(df, "variable1", "variable2", "group", "X-axis Label", "Y-axis Label", custom_colors, save_plot = FALSE)