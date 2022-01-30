#Figure 5: FISH image analysis
#July2021

library(tidyverse)

#read in table of nuclei counts per cluster per image
FISHtable = read.csv(file = "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/extractedData/FISH_data/FISHtestallprobes_clustercounts.csv")[ ,1:10]

plotDirectory = "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Figure5/"

#walk through each gene individually
#TNNT2
FISH_TNNT2 = FISHtable %>% filter(TNNT2>0)
TNNT2hist = ggplot(data=FISH_TNNT2, aes(x=TNNT2))+
  geom_histogram(fill = 'green',color = 'grey95',binwidth = 5, boundary=0.5)+
  geom_vline(xintercept=mean(FISH_TNNT2$TNNT2),color='darkgreen') +
  geom_text(aes(x=mean(TNNT2)+100, y=40, label=paste0('mean cluster size = ',round(mean(TNNT2),0))))+
  coord_cartesian(xlim = c(0, 400), ylim=c(0,55))+
  ggtitle(paste0(nrow(FISH_TNNT2)," TNNT2-high clusters from ",nrow(FISH_TNNT2 %>% select(Well,Image) %>% unique())," images"))+
  xlab("cluster size")+ylab("count")+
  theme_bw()
#save histogram of TNNT2-high cell clusters
pdf(paste(plotDirectory,"TNNT2hist.pdf",sep=""), width=6, height=4)
print(TNNT2hist)
dev.off()

#cell-based analysis
uncountFISH_TNNT2 = FISH_TNNT2 %>% 
  select(Well, Image, Cluster, TNNT2) %>% #remove other columns
  uncount() #cell-based
summary(uncountFISH_TNNT2$TNNT2)
#proportion of TNNT2-high cells found in clusters of >20
nrow(filter(uncountFISH_TNNT2, TNNT2>20))/nrow(uncountFISH_TNNT2)


#LUM
FISH_LUM = FISHtable %>% filter(LUM>0)
LUMhist = ggplot(data=FISH_LUM, aes(x=LUM))+
  geom_histogram(fill = 'red',color = 'grey95',binwidth = 5, boundary=0.5)+
  geom_vline(xintercept=mean(FISH_LUM$LUM),color='darkred') +
  geom_text(aes(x=mean(LUM)+100, y=40, label=paste0('mean cluster size = ',round(mean(LUM),0))))+
  coord_cartesian(xlim = c(0, 400), ylim=c(0,55))+
  ggtitle(paste0(nrow(FISH_LUM)," LUM-high clusters from ",nrow(FISH_LUM %>% select(Well,Image) %>% unique())," images"))+
  xlab("cluster size")+ylab("count")+
  theme_bw()
#save histogram of LUM-high cell clusters
pdf(paste(plotDirectory,"LUMhist.pdf",sep=""), width=6, height=4)
print(LUMhist)
dev.off()

#cell-based
uncountFISH_LUM = FISH_LUM %>% 
  select(Well, Image, Cluster, LUM) %>% #remove other columns
  uncount() #cell-based
summary(uncountFISH_LUM$LUM)
#proportion of LUM-high cells found in clusters of >20
nrow(filter(uncountFISH_LUM, LUM>20))/nrow(uncountFISH_LUM)


#EPCAM
FISH_EPCAM = FISHtable %>% filter(EPCAM>0)
EPCAMhist = ggplot(data=FISH_EPCAM, aes(x=EPCAM))+
  geom_histogram(fill = 'orange',color = 'grey95',binwidth = 5, boundary=0.5)+
  geom_vline(xintercept=mean(FISH_EPCAM$EPCAM),color='darkorange') +
  geom_text(aes(x=mean(EPCAM)+100, y=40, label=paste0('mean cluster size = ',round(mean(EPCAM),0))))+
  coord_cartesian(xlim = c(0, 400), ylim=c(0,55))+
  ggtitle(paste0(nrow(FISH_EPCAM)," EPCAM-high clusters from ",nrow(FISH_EPCAM %>% select(Well,Image) %>% unique())," images"))+
  xlab("cluster size")+ylab("count")+
  theme_bw()
#save histogram of EPCAM-high cell clusters
pdf(paste(plotDirectory,"EPCAMhist.pdf",sep=""), width=6, height=4)
print(EPCAMhist)
dev.off()

#cell-based
uncountFISH_EPCAM = FISH_EPCAM %>% 
  select(Well, Image, Cluster, EPCAM) %>% #remove other columns
  uncount() #cell-based
summary(uncountFISH_EPCAM$EPCAM)
#proportion of EPCAM-high cells found in clusters of >20
nrow(filter(uncountFISH_EPCAM, EPCAM>20))/nrow(uncountFISH_EPCAM)


#ISL1
FISH_ISL1 = FISHtable %>% filter(ISL1>0)
ISL1hist = ggplot(data=FISH_ISL1, aes(x=ISL1))+
  geom_histogram(fill = 'pink',color = 'grey95',binwidth = 5, boundary=0.5)+
  geom_vline(xintercept=mean(FISH_ISL1$ISL1),color='darkmagenta') +
  geom_text(aes(x=mean(ISL1)+100, y=40, label=paste0('mean cluster size = ',round(mean(ISL1),0))))+
  coord_cartesian(xlim = c(0, 400), ylim=c(0,55))+
  ggtitle(paste0(nrow(FISH_ISL1)," ISL1-high clusters from ",nrow(FISH_ISL1 %>% select(Well,Image) %>% unique())," images"))+
  xlab("cluster size")+ylab("count")+
  theme_bw()
#save histogram of ISL1-high cell clusters
pdf(paste(plotDirectory,"ISL1hist.pdf",sep=""), width=6, height=4)
print(ISL1hist)
dev.off()

#cell-based
uncountFISH_ISL1 = FISH_ISL1 %>% 
  select(Well, Image, Cluster, ISL1) %>% #remove other columns
  uncount() #cell-based
summary(uncountFISH_ISL1$ISL1)
#proportion of ISL1-high cells found in clusters of >20
nrow(filter(uncountFISH_ISL1, ISL1>20))/nrow(uncountFISH_ISL1)


#WT1
FISH_WT1 = FISHtable %>% filter(WT1>0)
WT1hist = ggplot(data=FISH_WT1, aes(x=WT1))+
  geom_histogram(fill = 'yellow',color = 'grey95',binwidth = 5, boundary=0.5)+
  geom_vline(xintercept=mean(FISH_WT1$WT1),color='tan') +
  geom_text(aes(x=mean(WT1)+100, y=40, label=paste0('mean cluster size = ',round(mean(WT1),0))))+
  coord_cartesian(xlim = c(0, 400), ylim=c(0,55))+
  ggtitle(paste0(nrow(FISH_WT1)," WT1-high clusters from ",nrow(FISH_WT1 %>% select(Well,Image) %>% unique())," images"))+
  xlab("cluster size")+ylab("count")+
  theme_bw()
#save histogram of WT1-high cell clusters
pdf(paste(plotDirectory,"WT1hist.pdf",sep=""), width=6, height=4)
print(WT1hist)
dev.off()

#cell-based
uncountFISH_WT1 = FISH_WT1 %>% 
  select(Well, Image, Cluster, WT1) %>% #remove other columns
  uncount() #cell-based
uncountWT1hist = ggplot(data=uncountFISH_WT1, aes(x=WT1))+
  geom_histogram(fill = 'yellow',color = 'grey95',binwidth = 5, boundary=0.5)+
  coord_cartesian(xlim = c(0, 400), ylim=c(0,450))+
  ggtitle(paste0(nrow(uncountFISH_WT1)," WT1-high cells in \n",nrow(FISH_WT1)," clusters from ",nrow(FISH_WT1 %>% select(Well,Image) %>% unique())," images"))+
  xlab("cluster size")+ylab("number of cells")+
  theme_bw()
summary(uncountFISH_WT1$WT1)
#proportion of WT1-high cells found in clusters of >20
nrow(filter(uncountFISH_WT1, WT1>20))/nrow(uncountFISH_WT1)

