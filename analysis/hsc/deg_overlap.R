# DEG overlap--------
## Venn
library(VennDiagram)

DI_I <- fread("G:/radiation/result/17_DEG_HSPC/DMSO-IR_IR.csv") %>% as.data.frame()
colnames(DI_I)[2:6] <- paste0("DI_I.",colnames(DI_I)[2:6])

I_C <- fread("G:/radiation/result/17_DEG_HSPC/IR_Ctrl.csv") %>% as.data.frame()
colnames(I_C)[2:6] <- paste0("I_C.",colnames(I_C)[2:6])

D_C <- fread("G:/radiation/result/17_DEG_HSPC/DMSO_Ctrl.csv") %>% as.data.frame()
colnames(D_C)[2:6] <- paste0("D_C.",colnames(D_C)[2:6])

m = merge(DI_I,I_C,by = "V1")
m = merge(m,D_C,by="V1")

Rps = m[grep("^Rp[sl]",m$V1),'V1']
m = m[!(m$V1 %in% Rps),]

mm = m[m$DI_I.p_val_adj < 0.05 & m$D_C.p_val_adj < 0.05 & m$I_C.p_val_adj < 0.05,]

I <- m[m$I_C.p_val_adj < 0.05,'V1']
D <- m[m$D_C.p_val_adj < 0.05,'V1']
DI <- m[m$DI_I.p_val_adj < 0.05,'V1']
DEGs <- c(list(I),list(D),list(DI))
names(DEGs) <- c('IR vs Ctrl','DMSO vs Ctrl','DMSO-IR vs IR')

p = venn.diagram(
  x = DEGs,
  filename=NULL,
  fill= c("DarkRed","OliveDrab","RoyalBlue"),
  width = 1000,
  height = 1000, 
)

grid.draw(p)

