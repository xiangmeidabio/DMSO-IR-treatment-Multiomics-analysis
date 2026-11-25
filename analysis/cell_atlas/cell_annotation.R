# Marker------------
library(ggplot2)
ss = 1
t1 = FeaturePlot(merge_inter,features =c('Cd34'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","DarkRed"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t2 = FeaturePlot(merge_inter,features =c('Kit'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","DarkRed"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t3 = FeaturePlot(merge_inter,features =c('Hlf'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","DarkRed"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)


t4 = FeaturePlot(merge_inter,features =c('Car1'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","Coral"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t5 = FeaturePlot(merge_inter,features =c('Pf4'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","Sienna"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t6 = FeaturePlot(merge_inter,features =c('Elane'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","OliveDrab"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t7 = FeaturePlot(merge_inter,features =c('Cebpe'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","RoyalBlue"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t8 = FeaturePlot(merge_inter,features =c('Ly86'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","ForestGreen"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)


t9 = FeaturePlot(merge_inter,features =c('Siglech'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","ForestGreen"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t10 = FeaturePlot(merge_inter,features =c('Cd74'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","DarkSlateGray"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t11 = FeaturePlot(merge_inter,features =c('Prss34'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","DarkCyan"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t12 = FeaturePlot(merge_inter,features =c('Prg2'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","CadetBlue"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t13 = FeaturePlot(merge_inter,features =c('Flt3'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","DarkOrchid"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t14 = FeaturePlot(merge_inter,features =c('Dntt'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","DarkOrchid"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)


t15 = FeaturePlot(merge_inter,features =c('Tcf7'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","MediumPurple"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t16 = FeaturePlot(merge_inter,features =c('Ccl5'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","Orchid"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t17 = FeaturePlot(merge_inter,features =c('Gata3'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","Plum"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t18 = FeaturePlot(merge_inter,features =c('Cd3e'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","MediumVioletRed"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t19 = FeaturePlot(merge_inter,features =c('Cd19'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","Red"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

t20 = FeaturePlot(merge_inter,features =c('Mki67'),pt.size =ss,reduction = 'tsne')&
  scale_color_gradientn(colors = colorRampPalette(c("grey","black"))(100)) &
  NoAxes() &
  theme(aspect.ratio = 1)

(t1 | t2 | t3 | t4 | t5) / (t6 | t7 | t8 | t9 | t10) / (t11 | t12 | t13 | t14 | t15) / (t16 | t17 | t18 | t19 | t20)

# sample correlation------------
a<-AverageExpression(merge_inter,group.by = "orig.ident",slot = "data")
a<-a[["RNA"]]
a_dense <- as.matrix(a) 
res <- cor(a_dense, method = "pearson")

rdwhbu <- colorRampPalette(c("navy", "white", "brown3"))

pheatmap(res,
         color = rdwhbu(100),
         fontsize=10)

pheatmap(res,
         display_numbers = TRUE,
         color = rdwhbu(100),
         fontsize=10)

# Ratio change------------
## Radar plot
library(ggradar)
library(tidyverse)

t<-table(merge_inter$celltype,merge_inter$orig.ident)

t <- prop.table(t,2)
t <- as.data.frame(t)

df <- spread(t,Var2,Freq)
summary(df)

df[,6:9] <- df[,2:5]/df[,2]

df <- t(df)
colnames(df) <- df[1,]
df <- df[6:9,]

df <- apply(df, 2, as.numeric)
df <- apply(df, 2, log2)
rownames(df) <- c("Ctrl","log2(IR/Ctrl)","log2(DMSO/Ctrl)","log2(DMSO-IR/Ctrl)")

df %>% as.data.frame() %>% rownames_to_column("group") -> df

df <- df[2:4,]
df$group <- factor(df$group,levels = c("log2(IR/Ctrl)","log2(DMSO/Ctrl)","log2(DMSO-IR/Ctrl)"))

colnames(df)[1] <- "group"

ggradar(
  df, 
  values.radar = c("-4", "0", "3"),
  grid.min = -4, grid.mid = 0, grid.max = 3,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = color4[2:4],
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom"
)

## Bar plot
library(tidyr)
library(ggplot2)
library(ddply)

m <- merge_inter@meta.data
m$number=1

ggplot(m,aes(orig.ident,number,fill=celltype))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 16)+
  scale_fill_manual(values = colors3)+
  RotatedAxis()

m <- ddply(m,'orig.ident',transform,percent = 1/sum(number)*100)
ggplot(m,aes(orig.ident,percent,fill=celltype))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 16)+
  scale_fill_manual(values = colors3)+
  RotatedAxis()

# create some dummy data

t<-table(merge_inter$celltype,merge_inter$orig.ident)
t <- prop.table(t,2)
tt <- as.data.frame(t)

# create a grouped bar plot with ggplot2
ggplot(tt, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Var1, nrow = 2) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Cell type", y = "Proportion", fill = "Cell type")+
  scale_fill_manual(values = colors3)


