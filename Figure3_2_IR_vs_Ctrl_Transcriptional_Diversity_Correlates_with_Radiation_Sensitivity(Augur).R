# 1. Environment setup: Load required R packages and data ----
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(ggrepel)

# Load the integrated Seurat object
merge_inter <- readRDS("G:/radiation/RDS/merge_inte_signiture.rds")

# 2. Data preparation: Annotate cell types and merge with Augur results ----
n_clusters <- length(unique(merge_inter@meta.data$integrated_snn_res.0.6))
celltype_annotation <- data.frame(
  ClusterID = 0:(n_clusters - 1),
  celltype = "unknown"
)

# Assign biological identities to each cluster
celltype_annotation[celltype_annotation$ClusterID %in% c(0), "celltype"] <- "NMP"
celltype_annotation[celltype_annotation$ClusterID %in% c(1), "celltype"] <- "NMP-Cycle"
celltype_annotation[celltype_annotation$ClusterID %in% c(2), "celltype"] <- "HSC"
celltype_annotation[celltype_annotation$ClusterID %in% c(3), "celltype"] <- "NP"
celltype_annotation[celltype_annotation$ClusterID %in% c(4), "celltype"] <- "MP"
celltype_annotation[celltype_annotation$ClusterID %in% c(5), "celltype"] <- "MkP"
celltype_annotation[celltype_annotation$ClusterID %in% c(6), "celltype"] <- "CLP"
celltype_annotation[celltype_annotation$ClusterID %in% c(7), "celltype"] <- "BP-3"
celltype_annotation[celltype_annotation$ClusterID %in% c(8), "celltype"] <- "BP-2"
celltype_annotation[celltype_annotation$ClusterID %in% c(9), "celltype"] <- "ILC2"
celltype_annotation[celltype_annotation$ClusterID %in% c(10), "celltype"] <- "BP-1"
celltype_annotation[celltype_annotation$ClusterID %in% c(11), "celltype"] <- "ILC1"
celltype_annotation[celltype_annotation$ClusterID %in% c(12), "celltype"] <- "HSC-Cycle"
celltype_annotation[celltype_annotation$ClusterID %in% c(13), "celltype"] <- "T Cell"
celltype_annotation[celltype_annotation$ClusterID %in% c(14), "celltype"] <- "EP"
celltype_annotation[celltype_annotation$ClusterID %in% c(15), "celltype"] <- "DC"
celltype_annotation[celltype_annotation$ClusterID %in% c(16), "celltype"] <- "ILCP"
celltype_annotation[celltype_annotation$ClusterID %in% c(17), "celltype"] <- "MPP"
celltype_annotation[celltype_annotation$ClusterID %in% c(18), "celltype"] <- "B Cell"
celltype_annotation[celltype_annotation$ClusterID %in% c(19), "celltype"] <- "PreE"

# Add cell type annotation to Seurat object metadata
merge_inter$celltype <- "unknown"
for (i in 1:nrow(celltype_annotation)) {
  idx <- which(merge_inter@meta.data$integrated_snn_res.0.6 == celltype_annotation$ClusterID[i])
  merge_inter@meta.data[idx, "celltype"] <- celltype_annotation$celltype[i]
}

# Verify all cells have been annotated
print("Cell type annotation summary:")
print(table(merge_inter$celltype))

# Load Augur results and prepare data frame for correlation analysis
augur_results <- readRDS("G:/radiation/result/6_AUGUR/Group/IR_Ctrl.rds")
auc_df <- as.data.frame(augur_results$AUC)
colnames(celltype_annotation) <- c("cell_type", "celltype")
tt <- merge(celltype_annotation, auc_df, by = "cell_type")
print("Final merged data frame 'tt' created:")
print(tt)

# 3. Define diversity metric calculation functions ----

# Function 1: Shannon entropy
calculate_shannon <- function(seurat_object, assay = "RNA", slot = "counts") {
  counts <- GetAssayData(object = seurat_object, assay = assay, slot = slot)
  n_umis <- Matrix::colSums(counts)
  props <- counts
  col_indices <- rep(1:ncol(counts), diff(counts@p))
  props@x <- props@x / n_umis[col_indices]
  props@x <- -props@x * log2(props@x)
  shannon_diversity <- Matrix::colSums(props)
  return(shannon_diversity)
}

# Function 2: Gini index
calculate_gini <- function(seurat_object, assay = "RNA", slot = "counts") {
  counts <- GetAssayData(object = seurat_object, assay = assay, slot = slot)
  gini_indices <- apply(counts, 2, function(cell_counts) {
    cell_counts <- cell_counts[cell_counts > 0]
    if (length(cell_counts) < 2) return(0)
    cell_counts <- sort(cell_counts)
    n <- length(cell_counts)
    gini <- (2 * sum(1:n * cell_counts) / (n * sum(cell_counts))) - (n + 1) / n
    return(gini)
  })
  return(gini_indices)
}

# 4. Compute diversity metrics and add to Seurat object ----
cat("Computing diversity metrics...\n")
merge_inter$shannon_diversity <- calculate_shannon(merge_inter)
merge_inter$gini_index <- calculate_gini(merge_inter)
cat("Diversity metrics computed successfully.\n")

# 5. Summarize median values of metrics in control group by cell type ----
meta_ctrl <- merge_inter@meta.data %>% filter(orig.ident == "Ctrl")
summary_stats <- meta_ctrl %>%
  group_by(celltype) %>%
  summarise(
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    median_shannon = median(shannon_diversity),
    median_gini = median(gini_index),
    .groups = "drop"
  )
print("Summary of median metrics in control group:")
print(summary_stats)

# 6. Perform correlation analysis between metrics and Augur AUC, then visualize ----

# Merge summary statistics with Augur AUC values
final_data <- merge(summary_stats, tt, by = "celltype")

# Define cell type order and corresponding colors
celltype_order <- c(
  "NMP", "NMP-Cycle", "HSC", "NP", "MP", "MkP", "CLP", "BP-3", "BP-2", "ILC2",
  "BP-1", "ILC1", "HSC-Cycle", "T Cell", "EP", "DC", "ILCP", "MPP", "B Cell", "PreE"
)

color_palette <- c(
  "DarkRed", "LightCoral", "Red", "Coral", "Sienna", "LightGreen", "OliveDrab",
  "RoyalBlue", "ForestGreen", "DarkSlateGray", "DarkCyan", "CadetBlue",
  "PaleTurquoise", "GoldEnrod", "DarkOrchid", "MediumPurple", "Orchid",
  "Plum", "MediumVioletRed", "pink"
)

# Ensure celltype is a factor with consistent levels
merge_inter$celltype <- factor(merge_inter$celltype, levels = celltype_order)
named_colors <- setNames(color_palette, celltype_order)

# Function to generate correlation plot
create_correlation_plot <- function(data, metric_col, auc_col = "auc", color_vector) {
  # Perform Spearman correlation test
  corr_test <- cor.test(data[[metric_col]], data[[auc_col]], method = "spearman")
  rho <- round(corr_test$estimate, 3)
  pval <- format.pval(corr_test$p.value, digits = 2, eps = 0.001)
  
  # Clean metric name for labeling
  clean_metric_name <- gsub("median_", "", metric_col)
  x_label <- paste("Median", clean_metric_name, "in Control Group")
  
  # Create plot
  p <- ggplot(data, aes(x = .data[[metric_col]], y = .data[[auc_col]])) +
    geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80", linewidth = 0.7) +
    geom_point(aes(color = celltype), size = 4, alpha = 0.9) +
    scale_color_manual(values = color_vector) +
    geom_text_repel(
      aes(label = celltype),
      size = 3.5,
      max.overlaps = 15,
      segment.color = "grey50",
      show.legend = FALSE
    ) +
    stat_cor(method = "spearman", label.x.npc = 0.05, label.y.npc = 0.98, size = 4.5) +
    labs(
      title = paste("Metric:", clean_metric_name),
      x = x_label,
      y = "Augur AUC (IR vs. Ctrl)"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  
  return(p)
}

# List of metrics to test
metrics_to_test <- c("median_nFeature", "median_nCount", "median_shannon", "median_gini")

# Initialize containers for results
plot_list <- list()
correlation_results <- data.frame()

# Loop over metrics: compute correlation and generate plots
for (metric in metrics_to_test) {
  corr_test <- cor.test(final_data[[metric]], final_data$auc, method = "spearman")
  correlation_results <- rbind(correlation_results, data.frame(
    Metric = metric,
    Correlation = corr_test$estimate,
    P_value = corr_test$p.value
  ))
  plot_list[[metric]] <- create_correlation_plot(
    data = final_data,
    metric_col = metric,
    color_vector = named_colors
  )
}

# 7. Output final resultÙÕù ¡c×   ƒ—ñû  ]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             òÍª,pÓ«©67^ä´ï+6o0…e7âºâ„MFõ%C­Zƒ°WocjÇÁrŒø'Þ:â!Hüh³°é3¢$ó"DIšU}sÀ¨AÝ°ì¡Ûõësj/œv51<þ‘ääú¦#š