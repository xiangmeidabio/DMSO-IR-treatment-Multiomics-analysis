## heatmap without RPS---------
library(Seurat)
library(ggplot2)
library(pheatmap)

HSC <- readRDS("G:/radiation/RDS/HSPC.lineage.rds")

HSC$bias = HSC$predicted.id
HSC@meta.data <- HSC@meta.data  %>%
  mutate(bias = case_when(
    bias %in% c("BP-1") ~ "Bas-bias",
    bias %in% c("NMP-1", "NMP-cycle") ~ "NM-bias",
    bias %in% c("CLP") ~ "Lym-bias",
    bias %in% c("MkP") ~ "Mk-bias",
    TRUE ~ bias))
HSC$bias = factor(HSC$bias,levels = c("Mk-bias","Bas-bias","NM-bias","Lym-bias"))

rgb_colors <- col2rgb(c("Coral","DarkCyan","LightGreen","DarkOrchid"))
alpha <- 0.4
new_colors <- rgb(rgb_colors[1,], rgb_colors[2,], rgb_colors[3,], alpha = alpha * 255, maxColorValue = 255)

DimPlot(HSC,group.by = "bias",
        cols = new_colors,
        pt.size = 3,
        reduction = "umap")+
  xlim(c(-5,5))

DI_I <- fread("G:/radiation/result/17_DEG_HSPC/DMSO-IR_IR.csv") %>% as.data.frame()
colnames(DI_I)[2:6] <- paste0("DI_I.",colnames(DI_I)[2:6])

I_C <- fread("G:/radiation/result/17_DEG_HSPC/IR_Ctrl.csv") %>% as.data.frame()
colnames(I_C)[2:6] <- paste0("I_C.",colnames(I_C)[2:6])

D_C <- fread("G:/radiation/result/17_DEG_HSPC/DMSO_Ctrl.csv") %>% as.data.frame()
colnames(D_C)[2:6] <- paste0("D_C.",colnames(D_C)[2:6])

m = merge(DI_I,I_C,by = "V1")
m = merge(m,D_C,by="V1")

mm = m[m$DI_I.p_val_adj < 0.05 & m$D_C.p_val_adj < 0.05 & m$I_C.p_val_adj < 0.05,]

Rps = mm[grep("^Rp[sl]",mm$V1),'V1']
mm = mm[!(mm$V1 %in% Rps),]

exp = AverageExpression(HSC,assays = "RNA",group.by = "orig.ident",features = mm$V1)
exp = exp$RNA
exp = as.matrix(exp)

p <- pheatmap(exp,scale = "row",cluster_cols = F,show_rownames = F)

t <- t(scale(t(exp)))

rdwhbu <- colorRampPalette(c("navy", "white", "brown3"))

p1 <- pheatmap(t,
               border_color=NA,
               color = rdwhbu(100),
               #scale = "row",
               #cutree_rows = N,
               cluster_row = T,
               cluster_col = F,
               show_rownames=T,
               show_colnames=T,
               clustering_distance_rows = "euclidean",
               clustering_method='complete',
               #breaks = breaksList,
               #annotation_row=annotation_row,
               #annotation_colors=ann_colors
)

## heatmap cluster-------
d = dist(t, method = 'euclidean')
tree = hclust(d, method = 'complete')

N=6
v = cutree(tree, N)[tree$order]
gaps = which((v[-1] - v[-length(v)]) != 0)

gene.cluster <- as.data.frame(v)
gene.cluster$gene <- rownames(gene.cluster)
table(gene.cluster$v)

annotation_row <- data.frame(Cut = gene.cluster$v
                             #gene = rownames(gene.cluster)
)

rownames(annotation_row) <- rownames(gene.cluster)

library(RColorBrewer)
c1 <- brewer.pal(N, "Set2")
names(c1) <- unique(annotation_row$Cut)

ann_colors = list(Cut = c1)

p <- pheatmap(t,
              border_color=NA,
              color = rdwhbu(100),
              #scale = "row",
              cutree_rows = N,
              cluster_row = T,
              cluster_col = F,
              show_rownames=F,
              show_colnames=T,
              clustering_distance_rows = "euclidean",clustering_method='complete',
              #breaks = breaksList,
              annotation_row=annotation_row,
              annotation_colors=ann_colors
)
p

