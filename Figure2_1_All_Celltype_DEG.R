## DEG---------
merge_inter$sample_celltype <- paste0(merge_inter$orig.ident,"_",merge_inter$celltype)
ct <- levels(merge_inter$celltype)
 
DefaultAssay(merge_inter) <- "RNA"
Idents(merge_inter) <- merge_inter$sample_celltype
 
for (t in list(c("IR","Ctrl"),c("DMSO","Ctrl"),c("DMSO-IR","IR"),c("DMSO-IR","DMSO"),c("DMSO-IR","Ctrl"))){
   t <- c("IR","Ctrl")
   for (i in ct ){
     i1 = paste0(t[1],"_",i)
     i2 = paste0(t[2],"_",i)
     print(paste0(i1," vs ",i2))
     tt <- FindMarkers(merge_inter,ident.1 = i1,ident.2 = i2)
 
     ii <- gsub("[/]","",i)
 
     fwrite(tt,paste0("G:/radiation/result/18_DEG_all_cluster/",t[1],"_",t[2],"_",ii,".csv"),row.names = T)
   }
 }

## Volcano each celltype---------
# merge DEG table
files = list.files('G:/radiation/result/18_DEG_all_cluster')
files <- files[grep("^IR_Ctrl", files)]
file_list = list()
for (f in files){
  deg = fread(paste0('G:/radiation/result/18_DEG_all_cluster/', f))
  c = sub("^IR_Ctrl_(.*)\\.csv$", "\\1", f)
  deg$cluster = c
  file_list = c(file_list,list(deg))
}
df <- do.call(rbind, file_list)

colnames(df)[1] = 'gene'

df <- df %>%
  mutate(label = case_when(
    p_val_adj < 0.05 & avg_log2FC > 0 ~ "Up",
    p_val_adj < 0.05 & avg_log2FC < 0 ~ "Down",
    TRUE ~ "Not Significant"
  ))
df$label = factor(df$label,levels = c('Up','Down',"Not Significant"))

df <- df  %>%
  mutate(cluster = case_when(
    cluster %in% c('ILC1PNKP') ~ "ILC1P/NKP",
    cluster %in% c("MP  pDC") ~ "MP / pDC",
    TRUE ~ cluster
  ))

df$cluster <- factor(df$cluster,levels = levels(merge_inter$celltype))

# top DEG
top10_list = list()
for (c in levels(df$cluster)){
  top10sig_up <- filter(df,cluster==c & p_val_adj < 0.05 & avg_log2FC >0) %>% distinct(gene,.keep_all = T) %>% top_n(3,avg_log2FC)
  top10sig_down <- filter(df,cluster==c  & p_val_adj < 0.05 & avg_log2FC <0) %>% distinct(gene,.keep_all = T) %>% top_n(3,-avg_log2FC)
  top10_list = c(top10_list,list(top10sig_up),list(top10sig_down))
}
top10_DEG <- do.call(rbind, top10_list)

df$size <- case_when(!(df$gene %in% top10_DEG$gene)~ 1,df$gene %in% top10_DEG$gene ~ 2)
dt <- filter(df,size==1)

# plot 
p <- ggplot()+
  geom_jitter(data = dt,aes(x = cluster, y = avg_log2FC, color = label), size = 0.85,width =0.4)+
  geom_jitter(data = top10_DEG, aes(x = cluster, y = avg_log2FC, color = label),size = 2,width =0.4)

# add cluster
dfcol<-data.frame(x=levels(df$cluster),y=0,label=levels(df$cluster))
p2 <- p + geom_tile(data = dfcol,aes(x=x,y=y),height=0.4,color = "black",fill = colors3[1:20],
                   alpha = 0.6,show.legend = F)

# add top gene name
library(ggplot2)
library(ggrepel)
p3 <- p2+
  geom_text_repel(data=top10_DEG,aes(x=cluster,y=avg_log2FC,label=gene),
                  # force = 1.2,
                  arrow = arrow(length = unit(0.008, "npc"),type = "open", ends = "last")
                  )

p4 <- p3 +scale_color_manual(name=NULL,values = c("brown","navy","grey"))

p5 <- p4+
  labs(x="Cluster",y="average logFC")+
  geom_text(data=dfcol,aes(x=x,y=y,label=label), size =3, color ="white")

p6 <- p5+
  theme_minimal()+
  theme(axis.title = element_text(size = 13,color = "black",face = "bold"),
  axis.line.y = element_line(color = "black", size = 1.2),axis.line.x = element_blank(),
  axis.text.x = element_blank(), panel.grid = element_blank(),legend.position = "top",
  legend.direction = "vertical",
  legend.justification = c(1,0),
  legend.text = element_text(size = 15)
  )

p6

## merged DEG matrix---------
t <- c("IR","Ctrl")
ct <- levels(merge_inter$celltype)
all.gene <- row.names(merge_inter)
 
result <- as.data.frame(matrix(nrow=length(all.gene),ncol = length(ct)+1,data = 0))
colnames(result) <- c("gene",ct)
result$gene <- all.gene
 
result_up <- result
result_down <- result

fc = 1
 
for (i in ct){
   print(i)
   i=ct[1]
   i1 = paste0(t[1],"_",i)
   i2 = paste0(t[2],"_",i)
   print(paste0(i1," vs ",i2))
 
   ii <- gsub("[/]","",i)
   tt <- fread(paste0("./result/18_DEG_all_cluster/",t[1],"_",t[2],"_",ii,".csv"))
 
 
   up <- tt[tt$p_val_adj < 0.05 & tt$avg_log2FC > fc,]
   up <- up$V1
 
   down <- tt[tt$p_val_adj < 0.05 & tt$avg_log2FC < -fc,]
   down <- down$V1
 
   result_up[result_up$gene %in% up,i] = 1
   result_down[result_down$gene %in% down,i] = 1
 
}
 
result_up$sum <- rowSums(result_up[2:ncol(result_up)])
result_down$sum <- rowSums(result_down[2:ncol(result_down)])
fwrite(result_up,"./result/19_HSPC_specific_DEG/UP.all.cluster_fc1.csv")
fwrite(result_down,"./result/19_HSPC_specific_DEG/DOWN.all.cluster_fc1.csv")
