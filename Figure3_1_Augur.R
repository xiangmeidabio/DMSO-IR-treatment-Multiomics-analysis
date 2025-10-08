# Augur
n=length(unique(merge_inter@meta.data$integrated_snn_res.0.6))
celltype=data.frame(ClusterID=0:(n-1),
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0),2]='NMP'
celltype[celltype$ClusterID %in% c(1),2]='NMP-Cycle'
celltype[celltype$ClusterID %in% c(2),2]='HSC'
celltype[celltype$ClusterID %in% c(3),2]='NP'
celltype[celltype$ClusterID %in% c(4),2]='MP' 
celltype[celltype$ClusterID %in% c(5),2]='MkP' 
celltype[celltype$ClusterID %in% c(6),2]='CLP'
celltype[celltype$ClusterID %in% c(7),2]='BP-3'
celltype[celltype$ClusterID %in% c(8),2]="BP-2"
celltype[celltype$ClusterID %in% c(9),2]="ILC2"
celltype[celltype$ClusterID %in% c(10),2]="BP-1"
celltype[celltype$ClusterID %in% c(11),2]="ILC1"
celltype[celltype$ClusterID %in% c(12),2]="HSC-Cycle"
celltype[celltype$ClusterID %in% c(13),2]="T Cell"
celltype[celltype$ClusterID %in% c(14),2]="EP"
celltype[celltype$ClusterID %in% c(15),2]='DC'
celltype[celltype$ClusterID %in% c(16),2]="ILCP"
celltype[celltype$ClusterID %in% c(17),2]="MPP"
celltype[celltype$ClusterID %in% c(18),2]="B Cell"
celltype[celltype$ClusterID %in% c(19),2]="PreE"

merge_inter@meta.data$celltype = "unknown"
for(i in 1:nrow(celltype)){
  merge_inter@meta.data[which(merge_inter@meta.data$integrated_snn_res.0.6 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(merge_inter@meta.data$celltype)

#IR vs Ctrl
i="G:/radiation/result/6_AUGUR/Group/IR_Ctrl"
augur <- readRDS(paste0(i,".rds"))
t <- as.data.frame(augur$AUC)

for(i in unique(t$cell_type)){
  auc <- t[t$cell_type==i,'auc']
  merge_inter@meta.data[which(merge_inter@meta.data$integrated_snn_res.0.6 == i),'auc'] <- auc
}

Idents(merge_inter) <- merge_inter@meta.data$orig.ident
Ctrl <- subset(merge_inter,ident="Ctrl")
IR <- subset(merge_inter,ident="IR")

colnames(celltype) <- c("cell_type","celltype")

tt <- merge(celltype,t,by="cell_type")
fwrite(tt,"G:/radiation/PDF/PDF_0810/11_augur_IR_Ctrl.csv")

tt <- tt[order(tt$auc),]
tt$celltype <- factor(tt$celltype,levels = tt$celltype)

colors3 <- c("DarkRed","LightCoral","Red",
             "Coral","Sienna",
             "LightGreen","OliveDrab","RoyalBlue","ForestGreen","DarkSlateGray",
             "DarkCyan","CadetBlue","PaleTurquoise","GoldEnrod",
             "DarkOrchid","MediumPurple","Orchid","Plum",
             "MediumVioletRed","pink")

merge_inter@meta.data$celltype <- factor(merge_inter$celltype,levels=mylevel)
names(colors3) <- levels(merge_inter@meta.data$celltype)

ggplot(data=tt,mapping=aes(x=celltype,y=auc,fill=celltype))+
  geom_bar(stat="identity",width = 0.7)+
  coord_flip()+
  theme_classic()+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values = colors3[as.character(tt$celltype)])+
  NoLegend()

######
library(ggrepel)
library(ggpubr)

c <- merge_inter@meta.data
c <- c[c$orig.ident == "Ctrl",]

medians <- aggregate(nFeature_RNA~celltype,c,median)

medians <- medians[order(-medians$nFeature_RNA),]
# medians$celltype <- reorder(medians$celltype,medians$nFeature_RNA,median)

augur <- readRDS("G:/radiation/result/6_AUGUR/Group/IR_Ctrl.rds")
t <- as.data.frame(augur$AUC)

celltype
colnames(celltype) <- c("cell_type","celltype")
tt <- merge(celltype,t,by="cell_type")

m <- merge(medians,tt,by="celltype")

m$celltype <- factor(m$celltype,levels = mylevel)
m <- m[order(m$celltype),]

ggscatter(m,x="nFeature_RNA",y="auc",
          size = 5,
          add="reg.line",
          conf.int = T,
          color="celltype",
          palette  = colors3[c(1:4,6:20)],
          add.params = list(color="black",fill="grey")
)+
  stat_cor(method = "spearman",
           label.x = 1500,
           label.y=0.95)+
  geom_text_repel(
    data = m,
    aes(label = celltype),
    size = 4,
    segment.color = "black", show.legend = FALSE,
    max.overlaps=100)+
  coord_cartesian(xlim = c(1200,6000))

