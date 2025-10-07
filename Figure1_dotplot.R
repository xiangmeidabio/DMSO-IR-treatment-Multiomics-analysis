## dotplot-----------
# reverse the order
merge_inter@meta.data$celltype <- factor(merge_inter$celltype,levels=mylevel[length(mylevel):1])
f <- c("Cd34","Kit","Hlf",
       "Mki67",
       # "Ly6a","Procr","Fgd5","Hoxb5", 
       "Itga2b","Gp1bb","Pf4","Zfpm1",
       "Hba-a2","Car1","Klf1","Car2","Gata1",
       "Prtn3","F13a1","Ly6c2", 'Lyz2', 
       "Gstm1","Elane","Gfi1", "Cebpe","Irf8","Csf1r","Ly86","Siglech",
       "Cd74","H2-Aa","Cst3",
       "Ms4a2","Cpa3","Prss34","Mcpt8","Ifitm1","Gata2","Csf1",
       "Prg2","Prg3",
       "Dntt","Flt3",
       "Tcf7",
       "Ccl5","Ncr1","Klrb1c","Klrc1","Xcl1","Il2rb",
       "Lztfl1","Il1rl1","Gata3",
       "Vpreb1","Fcrla","Ebf1",
       "Cd3e","Cd8a"
)

DotPlot(merge_inter,features = f,group.by = "celltype")+
  scale_color_gradientn(colors = c(rdwhbu(100)[1:45],rdwhbu(100)[55:100]))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))
# RotatedAxis()

merge_inter@meta.data$celltype <- factor(merge_inter@meta.data$celltype ,levels=mylevel)
