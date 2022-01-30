#Supplemental figure:J SD sibling analysis on melanoma cell data from Yogesh Goyal
#May2021

library(tidyverse)
library(ggplot2)
library(reshape2)

#pull data
extractDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/extractedData/Yogo_data/"
Yogosample <- read_tsv(file = paste0(extractDirectory,"cellID_twoexamples.tsv")) %>% select(-row)
Yogoclusters <- read_tsv(file = paste0(extractDirectory,"umapClusters_onlyBarcodedCells.tsv")) %>% select(-row)

plotDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Supplement/YogoJSD/"
#Jensen-Shannon divergence: use https://enterotype.embl.de/enterotypes.html
KLD <- function(x,y) sum(x * log2(x/y))
JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))

#convert to cellsperclustperSampleNum tibble form
#find baseline cellspercluster
allbaseline <- Yogoclusters %>% 
  group_by(seurat_clusters,sampleNum) %>%
  summarise(allbarcellsperclust = n())
sib1base <- allbaseline %>% filter(sampleNum == 'S1') #19 clusters

sib2base <- allbaseline %>% filter(sampleNum == 'S2')

addbarNum <- Yogosample %>% group_by(BC50StarcodeD8) %>%
  summarize(barCount=n()) %>%
  arrange(desc(barCount)) %>%
  mutate(barNum = paste0("B",row_number())) %>%
  transmute(BC50StarcodeD8, barNum)

Yogodata <- Yogosample %>%
  inner_join(addbarNum,by='BC50StarcodeD8')%>%
  select(barNum,cellID,sampleNum) %>%
  inner_join(Yogoclusters, by=c('cellID','sampleNum')) 

obscellsperclustYogo <- Yogodata %>%
  group_by(barNum,sampleNum,seurat_clusters) %>%
  summarize(barSamclustCount=n())

  
sib1obs_B1 <- obscellsperclustYogo %>% 
  filter(sampleNum == 'S1') %>% filter(barNum == 'B1') %>%
  full_join(sib1base,by = c('seurat_clusters','sampleNum')) %>%
  complete(seurat_clusters,sampleNum, fill = c(list(barNum = 'B1'),list(barSamclustCount= 0))) 

sib1obs_B2 <- obscellsperclustYogo %>% 
  filter(sampleNum == 'S1') %>% filter(barNum == 'B2') %>%
  full_join(sib1base,by = c('seurat_clusters','sampleNum')) %>%
  complete(seurat_clusters,sampleNum, fill = c(list(barNum = 'B2'),list(barSamclustCount= 0)))

sib2obs_B1 <- obscellsperclustYogo %>% 
  filter(sampleNum == 'S2') %>% filter(barNum == 'B1') %>%
  full_join(sib2base,by = c('seurat_clusters','sampleNum')) %>%
  complete(seurat_clusters,sampleNum, fill = c(list(barNum = 'B1'),list(barSamclustCount= 0))) 

sib2obs_B2 <- obscellsperclustYogo %>% 
  filter(sampleNum == 'S2') %>% filter(barNum == 'B2') %>%
  full_join(sib2base,by = c('seurat_clusters','sampleNum')) %>%
  complete(seurat_clusters,sampleNum, fill = c(list(barNum = 'B2'),list(barSamclustCount= 0))) 

  
addprop <- sib1obs_B1 %>%
  full_join(sib1obs_B2) %>%
  full_join(sib2obs_B1) %>%
  full_join(sib2obs_B2) %>% #19*4 = 76 rows
  mutate(propperbarclust = barSamclustCount/allbarcellsperclust) %>%
  ungroup() %>%
  group_by(barNum,sampleNum) %>%
  mutate(prop = propperbarclust/sum(propperbarclust))


####measure pairwise JS distance across unrelated and related siblings
#make a dataframe of all sib1 probabilities and a dataframe of all sib2 probabilities
#then run a for loop through each table?
propsib1 = addprop %>%
  filter(sampleNum == 'S1') %>%
  ungroup() %>%
  mutate(prob = prop+0.000001) %>% #add to avoid NaNs due to zeroes
  select(barNum,seurat_clusters,prob) %>%
  pivot_wider(id_cols = seurat_clusters,names_prefix = "prob1",names_from = barNum,values_from = prob)

propsib2 = addprop %>%
  filter(sampleNum == 'S2') %>%
  ungroup() %>%
  mutate(prob = prop+0.000001) %>% #add to avoid NaNs due to zeroes
  select(barNum,seurat_clusters,prob) %>%
  pivot_wider(id_cols = seurat_clusters,names_prefix = "prob2",names_from = barNum,values_from = prob)

JSDMatrix = matrix(, nrow = 2, ncol = 2)
#generates empty matrix that is 18x18 (for the 18 barcodes that are associated with at least 5 cells in both arms)
for (i in c(2:length(propsib1))) {
  for (j in c(2:length(propsib2))) {
    JSDMatrix[i-1,j-1] = JSD(propsib1[,i],propsib2[,j]) #start at 2 so have to index at i-1, j-1
  }
}

usenames = unique(addprop$barNum)
JSDMatrixplot = JSDMatrix
rownames(JSDMatrixplot) = sub("$", "A", usenames)#colnames(propsib1)[2:3]
colnames(JSDMatrixplot) = sub("$", "B", usenames)#colnames(propsib2)[2:3]
melted_JSDMatrix <- melt(JSDMatrixplot, na.rm = TRUE) #list of JSD values
melted_JSDMatrix$value <- round(melted_JSDMatrix$value,2) #round to 2 decimal points
#try to reverse order
#https://stackoverflow.com/questions/40725146/how-do-i-reverse-the-y-axis-in-ggplots-geom-tile

#pull for rectanges
diag_JSDMatrix <- melted_JSDMatrix[which(melted_JSDMatrix$Var2 == sub("A", "B", melted_JSDMatrix$Var1)),]
#plots_JSDMatrix <- na.omit(diag_JSDMatrix[which(melted_JSDMatrix$Var1 %in% c("B14A","B18A","B23A")),])

plot_JSD = ggplot(data = melted_JSDMatrix, aes(x=Var2, y=forcats::fct_rev(Var1), fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient(low = "darkblue", high = "white", 
                      limit = c(0,1), 
                      breaks = seq(0,1,0.2),
                      labels = seq(0,1,0.2),
                      space = "Lab", 
                      name="Jensen-Shannon distance") +
  theme_minimal()+ 
  #https://stackoverflow.com/questions/13258454/marking-specific-tiles-in-geom-tile-geom-raster
  geom_rect(data=diag_JSDMatrix, size=1, fill=NA, colour="black", #used to use #FFFF6D
            aes(xmin=as.integer(Var2) - 0.5, xmax=as.integer(Var2) + 0.5, ymin=as.integer(forcats::fct_rev(Var1)) - 0.5, ymax=as.integer(forcats::fct_rev(Var1)) + 0.5)) +
  #geom_rect(data=plots_JSDMatrix, size=1, fill=NA, colour="#ff068f",
  #          aes(xmin=as.integer(Var2) - 0.5, xmax=as.integer(Var2) + 0.5, ymin=as.integer(forcats::fct_rev(Var1)) - 0.5, ymax=as.integer(forcats::fct_rev(Var1)) + 0.5)) +
  # theme(axis.text.x = element_text(angle = 30, vjust = 1, 
  #                                  size = 10, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Var2, forcats::fct_rev(Var1), label = value), color = "black", size = 3) +
  theme_classic((base_size = 10)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(colour = 'black', size = 1.5),
    axis.text.x = element_text(angle = 90),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_line(colour = "black", size = 1.5),
    #legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
    text=element_text(family="Helvetica"))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  scale_y_discrete(position = "right")
#theme(legend.position = "none")

#print(plot_JSD)
#save 700x700
pdf(paste(plotDirectory,"Yogo_JSD.pdf", sep=""), width=7, height=7)
print(plot_JSD)
dev.off()


addprop_B1 = addprop %>% filter(barNum == 'B1')
addprop_B2 = addprop %>% filter(barNum == 'B2')

samplenum.labs <- c(paste0("B1A",": ",sum(sib1obs_B1$barSamclustCount)," cells"),
                    paste0("B1B",": ",sum(sib2obs_B1$barSamclustCount)," cells"))
names(samplenum.labs) <- c("S1", "S2")
JS_B1 = JSD(propsib1[,2],propsib2[,2])
plotnormcellsclustYogo_B1 <- ggplot(addprop_B1,aes(x=seurat_clusters,y=prop,group=sampleNum,fill=sampleNum)) +
  geom_bar(stat='identity')+
  facet_wrap(~sampleNum,ncol=1,labeller = labeller(sampleNum = samplenum.labs))+
  ylab(paste0("normalized proportion of B1 cells"))+xlab("Seurat clusters, res=0.5")+
  scale_fill_manual(breaks=c("S1","S2"),values = c("#FFB6DB","#ff068f")) +
  ggtitle(paste("B1A vs B1B\nJensen-Shannon distance = ",round(JS_B1,3),sep=""))+
  theme_bw()+
  theme(legend.position = "none")
##save probability density bargraphs for Yogo barcode 1
pdf(paste0(plotDirectory,"Yogo_B1.pdf"), width=4, height=4)
print(plotnormcellsclustYogo_B1)
dev.off()

samplenum.labs <- c(paste0("B2A",": ",sum(sib1obs_B2$barSamclustCount)," cells"),
                    paste0("B2B",": ",sum(sib2obs_B2$barSamclustCount)," cells"))
names(samplenum.labs) <- c("S1", "S2")
JS_B2 = JSD(propsib1[,3],propsib2[,3])
plotnormcellsclustYogo_B2 <- ggplot(addprop_B2,aes(x=seurat_clusters,y=prop,group=sampleNum,fill=sampleNum)) +
  geom_bar(stat='identity')+
  facet_wrap(~sampleNum,ncol=1,labeller = labeller(sampleNum = samplenum.labs))+
  ylab(paste0("normalized proportion of B2 cells"))+xlab("Seurat clusters, res=0.5")+
  scale_fill_manual(breaks=c("S1","S2"),values = c("#FFB6DB","#ff068f")) +
  ggtitle(paste("B2A vs B2B\nJensen-Shannon distance = ",round(JS_B2,3),sep=""))+
  theme_bw()+
  theme(legend.position = "none")
#save probability density bargraphs for Yogo barcode 2 
pdf(paste0(plotDirectory,"Yogo_B2.pdf"), width=4, height=4)
print(plotnormcellsclustYogo_B2)
dev.off()


