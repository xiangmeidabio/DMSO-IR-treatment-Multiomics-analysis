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


