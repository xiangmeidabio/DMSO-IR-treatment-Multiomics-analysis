#### Homer plot new up---------
b <- fread("G:\radiation\ATAC\6_diffpeak_filter\Homer_IRup_FC1\knownResults.txt") %>% as.data.frame()

b <- b[b$`q-value (Benjamini)` < 0.05,]
colnames(b) <- c("Motif Name","Consensus","P-value","LogP","q-value (Benjamini)",
                 "Target Sequences with Motif","Ratio of Target Sequences with Motif",
                 "Background Sequences with Motif","Ratio of Background Sequences with Motif")
b$class <- "down"
b <- separate(b,col = "Motif Name",sep = "[(]",into=c("Motif","motif_info"),remove = F)

b$Motif <- toupper(b$Motif)

b$Ratio_in_target = as.numeric(gsub("%","",b$`Ratio of Target Sequences with Motif`))
b$FC <- (as.numeric(gsub("%","",b$`Ratio of Target Sequences with Motif`)))  / (as.numeric(gsub("%","",b$`Ratio of Background Sequences with Motif`))) 

b = b[!duplicated(b$Motif),]

b$Rank = 1:nrow(b)

# b[b$FC >=10,'FC'] = 10

ggplot(b,aes(x=Rank,y=-LogP,color=FC))+
  geom_point()+
  # scale_color_gradientn(colors = rdwhbu(100)[c(35:1)],limits = c(1, 10))+
  scale_color_gradientn(colors = rdwhbu(100)[c(70:100)])+
  geom_text_repel(
    data = b[1:15,],
    aes(label = Motif,colors='black'),
    size = 3,
    segment.color = "black", show.legend = FALSE,
    max.overlaps=100)+
  theme_bw()+
  theme( legend.title = element_blank(),)+
  ylab('-log10(p-value)')+
  xlab('Rank')

#### Homer plot new down---------
b <- fread("G:\radiation\ATAC\6_diffpeak_filter\Homer_IRdown_FC1\knownResults.txt") %>% as.data.frame()

b <- b[b$`q-value (Benjamini)` < 0.05,]
colnames(b) <- c("Motif Name","Consensus","P-value","LogP","q-value (Benjamini)",
                 "Target Sequences with Motif","Ratio of Target Sequences with Motif",
                 "Background Sequences with Motif","Ratio of Background Sequences with Motif")
b$class <- "down"
b <- separate(b,col = "Motif Name",sep = "[(]",into=c("Motif","motif_info"),remove = F)

b$Motif <- toupper(b$Motif)

b$Ratio_in_target = as.numeric(gsub("%","",b$`Ratio of Target Sequences with Motif`))
b$FC <- (as.numeric(gsub("%","",b$`Ratio of Target Sequences with Motif`)))  / (as.numeric(gsub("%","",b$`Ratio of Background Sequences with Motif`))) 

b = b[!duplicated(b$Motif),]

b$Rank = 1:nrow(b)

ggplot(b,aes(x=Rank,y=-LogP,color=FC))+
  geom_point()+
  # scale_color_gradientn(colors = rdwhbu(100)[c(35:1)],limits = c(1, 10))+
  scale_color_gradientn(colors = rdwhbu(100)[c(35:1)])+
  geom_text_repel(
    data = b[1:15,],
    aes(label = Motif),
    size = 3,
    segment.color = "black", show.legend = FALSE,
    max.overlaps=100)+
  theme_bw()+
  theme( legend.title = element_blank(),)+
  ylab('-log10(p-value)')+
  xlab('Rank')

#### Homer plot new down，GATA and KLF---------
# new down 
b <- fread("G:\radiation\ATAC\6_diffpeak_filter\Homer_IRdown_FC1\knownResults.txt") %>% as.data.frame()

b <- b[b$`q-value (Benjamini)` < 0.05,]
colnames(b) <- c("Motif Name","Consensus","P-value","LogP","q-value (Benjamini)",
                 "Target Sequences with Motif","Ratio of Target Sequences with Motif",
                 "Background Sequences with Motif","Ratio of Background Sequences with Motif")
b$class <- "down"
b <- separate(b,col = "Motif Name",sep = "[(]",into=c("Motif","motif_info"),remove = F)

b$Motif <- toupper(b$Motif)

b$Ratio_in_target = as.numeric(gsub("%","",b$`Ratio of Target Sequences with Motif`))
b$FC <- (as.numeric(gsub("%","",b$`Ratio of Target Sequences with Motif`)))  / (as.numeric(gsub("%","",b$`Ratio of Background Sequences with Motif`))) 

b = b[!duplicated(b$Motif),]

b$Rank = 1:nrow(b)
b[b$FC >=3,'FC'] = 3

b$TF_family = 'Other'
b[grep("GATA",b$Motif),'TF_family'] = 'GATA'

b[grep("KLF",b$Motif),'TF_family'] = 'KLF'

ggplot(b,aes(x=Rank,y=-LogP,color=TF_family))+
  geom_point()+
  scale_color_manual(values = c("DarkRed","navy","Grey"))+
  geom_text_repel(
    data = b[1:60,][b$TF_family %in% c('GATA','KLF'),],
    aes(label = Motif),
    size = 3,
    segment.color = "black", show.legend = FALSE,
    max.overlaps=20)+
  theme_bw()+
  theme( legend.title = element_blank(),)+
  ylab('-log10(p-value)')+
  xlab('Rank')
