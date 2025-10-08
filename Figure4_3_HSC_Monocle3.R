# Monocle3
## HSC
library(monocle3)

expression_matrix1 <-GetAssayData(HSPC,assay = "RNA",slot = "data")
dim(expression_matrix1)
cell_metadata1 <- HSPC@meta.data
gene_annotation1 <- data.frame(gene_short_name=rownames(expression_matrix1))
rownames(gene_annotation1)<-rownames(expression_matrix1)
 
cds <- new_cell_data_set(expression_matrix1,
                          cell_metadata = cell_metadata1,
                          gene_metadata = gene_annotation1)
 
cds <- preprocess_cds(cds, num_dim = 50)
 
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype")
 
## import seurat embed
cds.embed <- cds@int_colData$reducedDims$UMAP
# int.embed <- Embeddings(HSPC, reduction = "tsne")
int.embed <- Embeddings(HSPC, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
colnames(cell_metadata1)
plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype")

cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) +
ggtitle("label by partitionID")
 
p1
p2
 
cds <- learn_graph(cds)
 
cds.embed <- cds@int_colData$reducedDims$UMAP %>% as.data.frame()
cc <- cds.embed[cds.embed$UMAP_1 < -2.5,] %>% rownames()
 
cds <- order_cells(cds,root_cells = cc)
 
saveRDS(cds,"G:/radiation/RDS/HSPC.Monocle3_2.rds")
p <- plot_cells(cds,
            color_cells_by = "bias",
            show_trajectory_graph=T,
            #trajectory_graph_color="grey11",
            #trajectory_graph_color="white",
            trajectory_graph_segment_size=1.5,
            label_cell_groups=F,
            #label_groups_by_cluster=T,
            #group_label_size=5,
            label_roots=F,
            label_leaves=F,
            label_branch_points=F,
            graph_label_size=5,
            cell_size=1)+
   scale_color_manual(values = new_colors)+
   xlim(c(-5,5))

cds <- readRDS("G:/radiation/RDS/HSPC.Monocle3_2.rds")

p <- plot_cells(cds,
                color_cells_by = "pseudotime",
                show_trajectory_graph=T,
                trajectory_graph_segment_size=1.5,
                label_cell_groups=F,
                label_roots=F,
                label_leaves=F,
                label_branch_points=F,
                graph_label_size=5,
                cell_size=1)+
  scale_fill_gradientn(colors = c("navy", "white", "brown3"))+
  xlim(c(-5,5))
p
