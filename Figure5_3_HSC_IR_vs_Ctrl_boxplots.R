## boxplot----------
library(ggpubr)

a <- fread("G:/radiation/result/17_DEG_HSPC/IR_Ctrl.csv") %>% as.data.frame()
down <- a[a$avg_log2FC < -1 & a$p_val_adj < 0.05,'V1']

merge_inter  <- AddModuleScore(merge_inter,features = list(down),name = "IR down Score")

t <- merge_inter@meta.data

color4 <- c("#171717A0","DarkRed","OliveDrab","RoyalBlue")

p <- ggboxplot(t,x="orig.ident",y="IR down Score1",
               fill ='orig.ident',
               outlier.shape=NA)+
  scale_fill_manual(values = color4)+
  theme_classic()+
  RotatedAxis()+
  facet_wrap(~celltype, nrow = 2)
p

my_comparisons <- list(c("IR", "Ctrl"), c("DMSO", "Ctrl"),c("DMSO-IR", "Ctrl"),
                       c("DMSO", "IR"), c("DMSO-IR", "DMSO"),c("DMSO-IR", "IR"))

p <- ggboxplot(t, x = "orig.ident", y = "IR down Score1",
               fill = 'orig.ident', outlier.shape = NA) +
  scale_fill_manual(values = color4) +
  theme_classic() +
  RotatedAxis() +
  facet_wrap(~celltype, nrow = 2) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE,
    label.y = 1.5,
    step.increase = 0.1       
  )

p

## boxplot----------
a <- fread("G:/radiation/result/17_DEG_HSPC/IR_Ctrl.csv") %>% as.data.frame()
down <- a[a$avg_log2FC > 1 & a$p_val_adj < 0.05,'V1']
merge_inter  <- AddModuleScore(merge_inter,features = list(down),name = "IR up Score")

t <- merge_inter@meta.data

p <- ggboxplot(t,x="celltype",y="IR.up.Score1",
               fill ='orig.ident',outlier.shape=NA)+
  scale_fill_manual(values = color4)+
  RotatedAxis()


p <- ggboxplot(t, x = "orig.ident", y = "IR up Score1",
               fill = 'orig.ident', outlier.shape = NA) +
  scale_fill_manual(values = color4) +
  theme_classic() +
  RotatedAxis() +
  facet_wrap(~celltype, nrow = 2) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE,
    label.y = 1.2,
    step.increase = 0.1    
  )

p
