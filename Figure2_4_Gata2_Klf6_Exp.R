### GATA2 Klf6 exp
t = FetchData(merge_inter,vars = c(paste0("Gata",1:3),'Klf6',"celltype","orig.ident"))
t = t[t$orig.ident %in% c("Ctrl","IR"),]
# t = t[t$celltype %in% c("HSC","MkP","PreE","BP-1","BP-2","BP-3"),]

p <- ggboxplot(t, x = "orig.ident", y = "Gata2",
               fill = "orig.ident",
               outlier.shape = NA) +
  facet_wrap(~ celltype, nrow = 2) +
  scale_fill_manual(values = color4) +
  stat_compare_means(
    method = "t.test",  
    label = "p.signif",      
    # comparisons = list(c("组1", "组2")),  
    hide.ns = TRUE,      
    label.y = 3.5
  )
p

t = FetchData(merge_inter,vars = c(paste0("Gata",1:3),'Klf6',"celltype","orig.ident"))
t = t[t$orig.ident %in% c("Ctrl","IR"),]
p <- ggboxplot(t, x = "orig.ident", y = "Klf6",
               fill = "orig.ident",
               outlier.shape = NA) +
  facet_wrap(~ celltype, nrow = 2) +
  scale_fill_manual(values = color4) +
  stat_compare_means(
    method = "t.test", 
    label = "p.signif",    
    hide.ns = TRUE,  
    label.y = 4
  )
p