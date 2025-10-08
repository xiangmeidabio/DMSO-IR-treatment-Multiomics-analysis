### barplot
library(data.table)

a = fread("G:/radiation/result/17_DEG_HSPC/Cluster/GO_select.csv")
tf <- a[,c(3,6,11)]
colnames(tf)[1:2] <- c("GO terms","P-value")
tf <- tf[as.numeric(tf$`P-value`)<=0.05,]
tf$`P-value` <- -log10(as.numeric(tf$`P-value`))
rownames(tf) <- tf$`GO terms`

tf <- tf[nrow(tf):1,]

names(c1) <- 1:6

pdf("G:/radiation/result/17_DEG_HSPC/Cluster/GO_select.pdf",width = 13,height = 18)

p <- barplot(as.numeric(tf$`P-value`),horiz=T,
             xlim=c(0,max(tf$`P-value`)+0.5),axes=T,
             col=c1[tf$Class],
             xlab ="-log10(p-value)",
             cex.axis=1.3,cex.lab=1.5,border = NA)

for (iii in 1:nrow(tf)){
  p+text(0,(1.2*iii-0.6),tf$`GO terms`[iii],cex=1.6,pos=4)
}

dev.off()