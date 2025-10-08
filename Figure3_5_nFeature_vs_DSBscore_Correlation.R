# =============================================================================
# 1. Load libraries and integrated Seurat object
# =============================================================================

library(Seurat)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggplot2)
library(ggpubr)

merge_inter <- readRDS("G:/radiation/RDS/merge_inte_signiture.rds")
DefaultAssay(merge_inter) <- "RNA"

# =============================================================================
# 2. Define DNA repair gene set using GO terms
# =============================================================================

go_ids <- c("GO:0006281", "GO:0006302")  # DNA repair, double-strand break repair
dna_repair_genes <- unlist(mapIds(
  org.Mm.eg.db,
  keys = go_ids,
  column = "SYMBOL",
  keytype = "GO",
  multiVals = "list"
))
dna_repair_genes <- dna_repair_genes[!is.na(dna_repair_genes)]
dna_repair_list <- list(DSB_Repair = dna_repair_genes)

# =============================================================================
# 3. Compute module score for DNA repair activity
# =============================================================================

merge_inter <- AddModuleScore(
  object = merge_inter,
  features = dna_repair_list,
  name = "DSB_Repair_Score"
)

# =============================================================================
# 4. Prepare metadata for both groups (Ctrl and DMSO)
# =============================================================================

mylevel <- c(
  "HSC","HSC-Cycle","MPP",
  "MkP","PreE",
  "NMP","NMP-Cycle","NP","MP","DC",
  "BP-1", "BP-2" ,"BP-3","EP",
  "CLP", "ILCP","ILC1","ILC2",
  "T Cell","B Cell"
)

colors3 <- c("DarkRed","LightCoral","Red","Coral","Sienna","LightGreen","OliveDrab",
             "RoyalBlue","ForestGreen","DarkSlateGray","DarkCyan","CadetBlue","PaleTurquoise",
             "GoldEnrod","DarkOrchid","MediumPurple","Orchid","Plum","MediumVioletRed","pink")

merge_inter$celltype <- factor(merge_inter$celltype, levels = mylevel)
named_colors <- setNames(colors3, mylevel)

ctrl_metadata <- subset(merge_inter@meta.data, orig.ident == "Ctrl")
dmsa_metadata <- subset(merge_inter@meta.data, orig.ident == "DMSO")

# =============================================================================
# 5. Plot correlation between transcriptional complexity and DSB repair score for each group
# =============================================================================

plot_group <- function(metadata, group_name) {
  metadata$celltype <- factor(metadata$celltype, levels = mylevel)
  
  p_corr <- ggplot(metadata, aes(x = nFeature_RNA, y = DSB_Repair_Score1)) +
    geom_point(aes(color = celltype), alpha = 0.6, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    stat_cor(method = "spearman", label.x.npc = 0.05, label.y.npc = 0.95, size = 4) +
    scale_color_manual(values = setNames(colors3, mylevel), limits = mylevel) +
    theme_bw() +
    labs(title = paste0("Baseline Transcriptional Complexity vs. DSB Repair Activity (", group_name, ")"),
         x = "Transcriptional Complexity (nFeature_RNA)",
         y = "DSB Repair Module Score")
  
  print(p_corr)
  
  ggsave(paste0("G:/radiation/", group_name, "_nFeature_vs_DSBscore_Correlation_ggplot.pdf"), plot = p_corr, width = 8, height = 6)
}

# Generate plots for both groups
plot_group(ctrl_metadata, "Ctrl")
plot_group(dmsa_metadata, "DMSO")
