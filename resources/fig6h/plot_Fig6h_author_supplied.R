# Figure 6h plotting code supplied by the author on 26 November 2026.
# Requires the original Seurat object merge_inter in the R session.
# Analysis and plotting parameters are preserved as supplied.
# Archived for provenance; this script was not rerun during the manuscript audit.

library(Seurat)
library(tidyverse)
library(ggpubr)

DefaultAssay(merge_inter) <- "RNA"
genes_of_interest <- c("Uqcrb", "Ptma", "Rps2", "Gapdh", "Impa2", "Prtn3", "Ramp1", "Snrpf", "Cenpw")
expr_data <- FetchData(merge_inter, vars = c(genes_of_interest, "orig.ident"), layer = "data")
colnames(expr_data)[colnames(expr_data) == "orig.ident"] <- "Group"
target_groups <- c("Ctrl", "IR", "DMSO", "DMSO-IR")
plot_data <- expr_data %>% filter(Group %in% target_groups)
plot_data$Group <- factor(plot_data$Group, levels = target_groups)
plot_data_long <- plot_data %>% pivot_longer(cols = all_of(genes_of_interest), names_to = "Gene", values_to = "Expression")
plot_data_long$Gene <- factor(plot_data_long$Gene, levels = genes_of_interest)
my_comparisons <- list(c("Ctrl", "IR"), c("Ctrl", "DMSO"), c("Ctrl", "DMSO-IR"), c("IR", "DMSO-IR"))
group_colors <- c("Ctrl" = "#B61C28", "IR" = "#A7B2C3", "DMSO" = "#F4DFE1", "DMSO-IR" = "#D68790")
p <- ggplot(plot_data_long, aes(x = Group, y = Expression, fill = Group)) +
 geom_violin(trim = FALSE, alpha = 0.6, scale = "width") +
 geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8, color = "black") +
 stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", vjust = 0.5) +
 facet_wrap(~Gene, scales = "free_y", ncol = 5) +
 scale_fill_manual(values = group_colors) +
 theme_bw() +
 labs(y = "Normalized Expression", x = "", title = "Expression of 9 Genes in the Integrated Dataset") +
 theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"), axis.text.y = element_text(size = 10, color = "black"), strip.text = element_text(size = 10, face = "bold.italic"), plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), legend.position = "none")
print(p)
ggsave("9_Genes_Expression_All_Cells_Clean.pdf", plot = p, width = 12, height = 10)
