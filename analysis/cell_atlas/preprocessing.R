# R------------
library(Seurat)
library(patchwork)
library(dplyr)
library(export)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggsci)
library(clustree)
library(monocle)
library(monocle3)
library(data.table)
library(readxl)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(plyr)
library(SCENIC)
library(ggrepel)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)
library(scales)

# RDS------------
merge_inter <- readRDS("G:/radiation/RDS/merge_inte_signiture.rds")

# Celltype------------
n=length(unique(merge_inter@meta.data$integrated_snn_res.0.6))
celltype=data.frame(ClusterID=0:(n-1),
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0),2]='NMP'
celltype[celltype$ClusterID %in% c(1),2]='NMP-cycle'
celltype[celltype$ClusterID %in% c(2),2]='HSC'
celltype[celltype$ClusterID %in% c(3),2]='NP'
celltype[celltype$ClusterID %in% c(4),2]='MP' 
celltype[celltype$ClusterID %in% c(5),2]='MkP' 
celltype[celltype$ClusterID %in% c(6),2]='CLP'
celltype[celltype$ClusterID %in% c(7),2]='BP-3'
celltype[celltype$ClusterID %in% c(8),2]="BP-2"
celltype[celltype$ClusterID %in% c(9),2]="ILC2"
celltype[celltype$ClusterID %in% c(10),2]="BP-1"
celltype[celltype$ClusterID %in% c(11),2]="ILC1"
celltype[celltype$ClusterID %in% c(12),2]="HSC-Cycle"
celltype[celltype$ClusterID %in% c(13),2]="T Cell"
celltype[celltype$ClusterID %in% c(14),2]="EP"
celltype[celltype$ClusterID %in% c(15),2]='DC'
celltype[celltype$ClusterID %in% c(16),2]="ILCP"
celltype[celltype$ClusterID %in% c(17),2]="MPP"
celltype[celltype$ClusterID %in% c(18),2]="B Cell"
celltype[celltype$ClusterID %in% c(19),2]="PreE"

# celltype
merge_inter@meta.data$celltype = "unknown"
for(i in 1:nrow(celltype)){
  merge_inter@meta.data[which(merge_inter@meta.data$integrated_snn_res.0.6 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(merge_inter@meta.data$celltype)

Idents(merge_inter) <- merge_inter$celltype
unique(merge_inter$celltype)
mylevel <- c(
  "HSC","HSC-Cycle","MPP",
  "MkP","PreE",
  "NMP","NMP-cycle","NP","MP","DC",
  "BP-1", "BP-2" ,"BP-3","EP",
  "CLP", "ILCP","ILC1","ILC2",
  "T Cell","B Cell"
)
merge_inter@meta.data$celltype <- factor(merge_inter$celltype,levels=mylevel)

# Dimplot------------
## TSNE
rgb_colors <- col2rgb(colors3)
alpha <- 0.3
new_colors <- rgb(rgb_colors[1,], rgb_colors[2,], rgb_colors[3,], alpha = alpha * 255, maxColorValue = 255)

DimPlot(merge_inter,group.by = 'celltype',
        reduction = "tsne",
        pt.size = 1,
        cols = new_colors)

DimPlot(merge_inter,group.by = 'celltype',
        reduction = "tsne",
        pt.size = 1,
        split.by = "orig.ident",
        cols = new_colors)

## UMAP
DimPlot(merge_inter,group.by = 'celltype',
        reduction = "umap",
        pt.size = 1,
        #label = T,label.size = 5,
        cols = new_colors)

