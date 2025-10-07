## coserved UP DEG-------------------
library(org.Mm.eg.db)
result_up <- fread("./result/19_HSPC_specific_DEG/UP.all.cluster.csv") %>% as.data.frame()
cDEG <- result_up[result_up$sum >= 14,'gene']

ego <- enrichGO(
  gene          = cDEG,
  keyType       = "SYMBOL",
  OrgDb         =  org.Mm.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  #readable      = TRUE
)
t <- ego@result
t <- t[as.numeric(t$p.adjust )<=0.05,]
fwrite(t,"G:/radiation/result/result/19_HSPC_specific_DEG/2.Conserved_DEG.GO.csv")

p <- dotplot(ego)
p

## coserved Down DEG-------------------
library(org.Mm.eg.db)
result_down <- fread("G:/radiation/result/result/19_HSPC_specific_DEG/DOWN.all.cluster.csv") %>% as.data.frame()
cDEG <- result_down[result_down$sum >= 14,'gene']
ego <- enrichGO(
  gene          = cDEG,
  keyType       = "SYMBOL",
  OrgDb         =  org.Mm.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  #readable      = TRUE
)
t <- ego@result
t <- t[as.numeric(t$p.adjust )<=0.05,]
fwrite(t,"G:/radiation/result/result/19_HSPC_specific_DEG/2.Conserved_down.GO.csv")

p <- dotplot(ego)
p


