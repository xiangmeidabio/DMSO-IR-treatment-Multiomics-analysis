# =============================================================================
# Additional Analysis: Correlation between transcriptional diversity and 
# chemical perturbation sensitivity (DMSO vs. Ctrl)
# Note: Transcriptional diversity metrics were precomputed in the Seurat object.
# =============================================================================

# <<< Key Change 1: Load Augur results for DMSO vs. Ctrl comparison >>>
augur_dmso <- readRDS("G:/radiation/result/6_AUGUR/Group/DMSO_Ctrl.rds")
t_dmso <- as.data.frame(augur_dmso$AUC)
colnames(celltype_final) <- c("cell_type", "celltype")
tt_dmso <- merge(celltype_final, t_dmso, by = "cell_type")
cat("Successfully generated tt_dmso data frame for DMSO vs. Ctrl comparison:\n")
print(tt_dmso)

# <<< Key Change 2: Merge diversity summary with DMSO-specific Augur AUC >>>
final_data_dmso <- merge(summary_stats, tt_dmso, by = "celltype")

# <<< Key Change 3: Initialize containers for DMSO analysis results >>>
plot_list_dmso <- list()
correlation_results_dmso <- data.frame()

# Perform correlation analysis for each diversity metric
for (metric in metrics_to_test) {
  corr_test <- cor.test(final_data_dmso[[metric]], final_data_dmso$auc, method = "spearman")
  correlation_results_dmso <- rbind(correlation_results_dmso, data.frame(
    Metric = metric,
    Correlation = corr_test$estimate,
    P_value = corr_test$p.value
  ))
  
  # Generate plot with updated y-axis label
  p <- create_correlation_plot(
    data = final_data_dmso,
    metric_col = metric,
    color_vector = named_colors
  ) +
    labs(y = "Augur AUC (DMSO vs. Ctrl)")
  
  plot_list_dmso[[metric]] <- p
}

# =============================================================================
# Output results for DMSO vs. Ctrl analysis
# =============================================================================

# Display all correlation plots
for (metric in metrics_to_test) {
  print(plot_list_dmso[[metric]])
}

# Format and display correlation statistics
correlation_results_dmso <- correlation_results_dmso %>%
  mutate(Abs_Correlation = abs(Correlation)) %>%
  arrange(desc(Abs_Correlation))

cat("=================================================================\n")
cat("Correlation between cellular diversity metrics and Augur AUC\n")
cat("(DMSO vs. Ctrl comparison, Spearman correlation):\n")
print(correlation_results_dmso)
cat("=================================================================\n")

