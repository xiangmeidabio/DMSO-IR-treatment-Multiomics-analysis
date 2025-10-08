## HSC lineage bias Ratio
library(plyr)

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

m <- HSC@meta.data
m$number=1

bias_color = c("Coral","DarkCyan","LightGreen","DarkOrchid")

p1 <- ggplot(m,aes(orig.ident,number,fill=bias))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 16)+
  scale_fill_manual(values = bias_color)+
  RotatedAxis()
p1

m <- ddply(m,'orig.ident',transform,percent = 1/sum(number)*100)
p2 <- ggplot(m,aes(orig.ident,percent,fill=bias))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 16)+
  scale_fill_manual(values = bias_color)+
  RotatedAxis()
p2
