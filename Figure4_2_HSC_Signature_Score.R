# HSC Signature score

## stem score
#Nature 20 Single-cell lineage tracing unveils a role for TCF15 in haematopoiesis
library(readxl)
a <- read_excel("G:/radiation/public_data/Nature_20_TCF15_ext_FIG5_table_signature.xlsx")
a <- as.data.frame(a)

HSC_score <- list()
for (i in 1:length(a)){
   t <-  a[,i]
   t <- t[!is.na(t)]
   HSC_score <- c(HSC_score,list(t))
}
names(HSC_score) <- colnames(a)
 
DefaultAssay(merge_inter) <- "RNA"
merge_inter <- AddModuleScore(merge_inter,features = HSC_score,name = colnames(a))

f <-c("Stemscore_Giladi_et_al6")
FeaturePlot(HSC,features = f,
            reduction = "umap",
            pt.size = 3)&
  scale_color_gradientn(colors = rdwhbu(100)) &
  xlim(c(-5,5))

## Box Plot
library(ggpubr)
library(rstatix)

t <- HSC@meta.data
stat_wilcox <- t_test(group_by(t, bias), Stemscore_Giladi_et_al6 ~orig.ident,paired = F)  
stat_wilcox <- add_significance(stat_wilcox, 'p')  
stat_wilcox.test <-  add_xy_position(stat_wilcox, x = 'bias')

color4 <- c("#171717A0","DarkRed","OliveDrab","RoyalBlue")

ggboxplot(t, x = "bias", y = "Stemscore_Giladi_et_al6",
          fill = "orig.ident",
          palette = color4,
          outlier.shape = NA) +
  scale_fill_manual(values = color4) +
  RotatedAxis() +
  stat_pvalue_manual(stat_wilcox.test, 
                     #label = 'p', 
                     label = 'p.signif',
                     tip.length = 0.005,
                     hide.ns=T)

