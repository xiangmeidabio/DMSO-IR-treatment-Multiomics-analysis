#### 1.DEG HSC lym bias---------------
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

m <- HSC@meta.data
mm <- m[m$bias %in% c("Mk-bias","Lym-bias"),]
cc <- rownames(mm)
cds_sub <-  cds[,cc]

track_genes <- fread("G:/radiation/result/16_monocle3/1.track_genes_new_CLP.csv")

g <- subset(track_genes, q_value < 0.05 & morans_I > 0.25)
genes <- g$gene_short_name

pt.matrix <- exprs(cds_sub)[match(genes,rownames(rowData(cds_sub))),order(pseudotime(cds_sub))]

ttt <- as.data.frame(pt.matrix)

p <- pseudotime(cds_sub)
p <- p[order(p)]

cell <- mm[names(p),]

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;

l <- cell$bias
names(l) <- rownames(cell)

l <- factor(l,levels = c("Mk-bias","Lym-bias"))
l_color <- c("Coral","DarkOrchid")
names(l_color) <- c("Mk-bias","Lym-bias")

s <- cell$orig.ident
names(s) <- rownames(cell)
s <- factor(s,levels = c("Ctrl","IR","DMSO","DMSO-IR"))
s_color <- color4
names(s_color) <- c("Ctrl","IR","DMSO","DMSO-IR")

col_fun <- colorRamp2(
  c(0, 5, 10), 
  c("navy", "white", "brown3")
)

column_ha <- HeatmapAnnotation(
  lineage = l,
  sample = s,
  pseudotime = p,
  col=list(lineage = l_color,
           sample = s_color,
           pseudotime= col_fun))

htkm <- Heatmap(
  pt.matrix,
  name = "z-score",
  col = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize =8),
  km = 2,
  row_title_rot  = 0,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  top_annotation = column_ha
)
htkm

#### 1.1.get gene-----------
d <- dist(pt.matrix,method = "euclidean")
cluster2 <- hclust(d, method = "complete")
plot(cluster2,hang = -1,cex=0.6,axes=FALSE,ann=FALSE)
cut <- cutree(cluster2,2)
cut

gene1<- names(cut)[cut==1]
gene2<- names(cut)[cut==2]

#### 1.2.boxplot gene1----------
DefaultAssay(HSC) <- "RNA"
HSC2 <- subset(HSC,cells=cc)

tt <- AverageExpression(HSC2,features =gene1,group.by = c("orig.ident"),
                        slot = "data")
t_RNA <- as.data.frame(tt[["RNA"]])

t_RNA <- t_RNA[,1:2]
t_RNA$gene <- rownames(t_RNA)
t_long <- reshape2::melt(t_RNA, id.vars = c("gene"),
             value.name = "Exp")
colnames(t_long)[2] <- "Group" 

test = 'wilcox'

p <- ggpaired(t_long, 
              x = "Group", y = "Exp",
              color = "Group",  # 按基因着色
              line.color = "gray",
              line.size = 0.4,
              palette = "npg") +
  geom_text_repel(
    data = subset(t_long, Group == "Ctrl"),
    aes(label = gene), 
    #position = position_jitterdodge(),
    max.overlaps = 20)+
  scale_color_manual(values = color4)+
  stat_compare_means(method = test,paired = TRUE)
p

#### 1.3.boxplot gene2--------------
DefaultAssay(HSC) <- "RNA"
HSC2 <- subset(HSC,cells=cc)

tt <- AverageExpression(HSC2,features =gene2,group.by = c("orig.ident"),
                        slot = "data")
t_RNA <- as.data.frame(tt[["RNA"]])

t_RNA <- t_RNA[,1:2]
t_RNA$gene <- rownames(t_RNA)
t_long <- reshape2::melt(t_RNA, id.vars = c("gene"),
             value.name = "Exp")
colnames(t_long)[2] <- "Group" 


p <- ggpaired(t_long, 
              x = "Group", y = "Exp",
              color = "Group",  # 按基因着色
              line.color = "gray",
              line.size = 0.4,
              palette = "npg") +
  geom_text_repel(
    data = subset(t_long, Group == "Ctrl"),
    aes(label = gene), 
    #position = position_jitterdodge(),
    max.overlaps = 20)+
  scale_color_manual(values = color4)+
  stat_compare_means(method = test,paired = TRUE)
p

#### 2.DEG HSC mye lineage-------------
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

m <- HSC@meta.data
mm <- m[m$bias != "Lym-bias",]
cc <- rownames(mm)

cds_sub <-  cds[,cc]

track_genes <- fread("G:/radiation/result/16_monocle3/1.track_genes_new_myeloid.csv")

genes <- subset(track_genes, q_value < 0.05 & morans_I > 0.25)
genes <- genes$gene_short_name

pt.matrix <- exprs(cds_sub)[match(genes,rownames(rowData(cds_sub))),order(pseudotime(cds_sub))]

ttt <- as.data.frame(pt.matrix)

p <- pseudotime(cds_sub)
p <- p[order(p)]

cell <- mm[names(p),]

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;
#K means with 6 groups

l <- cell$bias
# l <- cell$integrated_snn_res.0.3
names(l) <- rownames(cell)

l <- factor(l,levels = c("Mk-bias","NM-bias","Bas-bias"))
l_color <- c("Coral","LightGreen","DarkCyan")
names(l_color) <- c("Mk-bias","NM-bias","Bas-bias")

s <- cell$orig.ident
names(s) <- rownames(cell)
s <- factor(s,levels = c("Ctrl","IR","DMSO","DMSO-IR"))
s_color <- color4
names(s_color) <- c("Ctrl","IR","DMSO","DMSO-IR")

col_fun <- colorRamp2(
  c(0, 5, 10), 
  c("navy", "white", "brown3")
)

column_ha <- HeatmapAnnotation(
  lineage = l,
  sample = s,
  pseudotime = p,
  col=list(lineage = l_color,
           sample = s_color,
           pseudotime= col_fun))

htkm <- Heatmap(
  pt.matrix,
  name = "z-score",
  col = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 6),
  km = 2,
  row_title_rot  = 0,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  top_annotation = column_ha
)
htkm

#### 2.1.get gene-----------
d <- dist(pt.matrix,method = "euclidean")
cluster2 <- hclust(d, method = "complete")
plot(cluster2,hang = -1,cex=0.6,axes=FALSE,ann=FALSE)
cut <- cutree(cluster2,2)
cut

gene1<- names(cut)[cut==1]
gene2<- names(cut)[cut==2]

#### 2.2.boxplot gene1-----------
DefaultAssay(HSC) <- "RNA"
HSC2 <- subset(HSC,cells=cc)

tt <- AverageExpression(HSC2, features =gene1, group.by = c("orig.ident"),
                        slot = "data")
t_RNA <- as.data.frame(tt[["RNA"]])

t_RNA <- t_RNA[,1:2]
t_RNA$gene <- rownames(t_RNA)
t_long <- reshape2::melt(t_RNA, id.vars = c("gene"),
             value.name = "Exp")
colnames(t_long)[2] <- "Group" 


p <- ggpaired(t_long, 
              x = "Group", y = "Exp",
              color = "Group",  # 按基因着色
              line.color = "gray",
              line.size = 0.4,
              palette = "npg") +
  geom_text_repel(
    data = subset(t_long, Group == "Ctrl"),
    aes(label = gene), 
    #position = position_jitterdodge(),
    max.overlaps = 20)+
  scale_color_manual(values = color4)+
  stat_compare_means(method = test,paired = TRUE)
p

#### 2.3.boxplot gene2-----------
DefaultAssay(HSC) <- "RNA"
HSC2 <- subset(HSC,cells=cc)

tt <- AverageExpression(HSC2,features =gene2,group.by = c("orig.ident"),
                        slot = "data")
t_RNA <- as.data.frame(tt[["RNA"]])

t_RNA <- t_RNA[,1:2]
t_RNA$gene <- rownames(t_RNA)
t_long <- reshape2::melt(t_RNA, id.vars = c("gene"),
             value.name = "Exp")
colnames(t_long)[2] <- "Group" 


p <- ggpaired(t_long, 
              x = "Group", y = "Exp",
              color = "Group",  # 按基因着色
              line.color = "gray",
              line.size = 0.4,
              palette = "npg") +
  geom_text_repel(
    data = subset(t_long, Group == "Ctrl"),
    aes(label = gene), 
    #position = position_jitterdodge(),
    max.overlaps = 20)+
  scale_color_manual(values = color4)+
  stat_compare_means(method = test,paired = TRUE)
p

