#Supplemental figure: generating "probability distributions"
#July2021

library(tidyverse)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(grid)
library(DescTools)

#pull Seurat info
SeuratDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/extractedData/createSeuratobj/"

#no clusters
plusneg_scTransform.integrated <- readRDS(file = paste0(SeuratDirectory,'plusneg_scTransform_50pcs_5000var_filter.rds'))

#clusters res=0.5
plusneg_scTransformclusters_res0.5.integrated = readRDS(file = paste0(SeuratDirectory,'plusneg_scTransform_50pcs_5000var_res0.5_filter.rds'))

###pull out log normalized counts, UMAP coordinates, and Seurat cluster info out of Seurat
plusneglogNormalizedCounts = plusneg_scTransform.integrated[['SCT']]@data #normalized log counts matrix (for counts only, replace "data" by "counts")
plusnegnormalizedCounts = plusneg_scTransform.integrated[['SCT']]@counts #normalized counts matrix, i.e. sequencing depth corrected counts
plusnegcells_count = colnames(plusnegnormalizedCounts) 
plusnegcells_count = sub("-1", "", plusnegcells_count)
plusnegcells_count_cellID = sub(".*_", "", plusnegcells_count)
plusnegcells_count_sampleNum = sub("_.*","",plusnegcells_count)
plusneglogNormalizedCounts = as_tibble(as.data.frame((t(as.matrix(plusneglogNormalizedCounts)))))
plusneglogNormalizedCounts = plusneglogNormalizedCounts %>% mutate(cellID = plusnegcells_count_cellID,
                                                                   SampleNum = plusnegcells_count_sampleNum,
                                                                   cellIDSN = plusnegcells_count)

plusnegumapCoordinates = (plusneg_scTransform.integrated[['umap']])@cell.embeddings #UMAP coordinates
plusnegcells_UMAP = rownames(plusnegumapCoordinates) #CellIds with Sample number as prefix
plusnegcells_UMAP = sub("-1", "", plusnegcells_UMAP)
plusnegcells_UMAP_cellID = sub(".*_", "", plusnegcells_UMAP)
plusnegcells_UMAP_sampleNum = sub("_.*","",plusnegcells_UMAP)
plusnegumapCoordinates = as_tibble(plusnegumapCoordinates) #UMAP coordinates in tibble
plusnegumapCoordinates = plusnegumapCoordinates %>% mutate(cellID = plusnegcells_UMAP_cellID,
                                                           SampleNum = plusnegcells_UMAP_sampleNum,
                                                           cellIDSN = plusnegcells_UMAP)

plusnegcells_clusters <- plusneg_scTransformclusters_res0.5.integrated@active.ident
plusnegcells_clusters <- as.data.frame(plusnegcells_clusters)
cluster_rowname <- rownames(plusnegcells_clusters)
cluster_cellID = sub(".*_", "", cluster_rowname)
cluster_cellID = sub("-1", "", cluster_cellID)
cluster_SN = sub("_.*","",cluster_rowname)
plusnegcells_clusters = as_tibble(plusnegcells_clusters)
plusnegcells_clusters <- plusnegcells_clusters %>% mutate(cellID = cluster_cellID,
                                                          SampleNum = cluster_SN)

plusneglogNormalizedCounts <- inner_join(plusneglogNormalizedCounts, plusnegcells_clusters, by = c("cellID", "SampleNum")) #17599
plusnegumapCoordinates <- inner_join(plusnegumapCoordinates, plusnegcells_clusters, by = c("cellID", "SampleNum")) #17599

#pull barcode info
barcodeDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/extractedData/createbarcodeobj/"
corelinCountTooverlaps30d6 = readRDS(file = paste0(barcodeDirectory,'corelinCountTooverlaps30d6.rds'))

###merging cellIDs and Barcodes with other gene data
plusneglogNormalizedCountsSubsetWBarcodes = inner_join(plusneglogNormalizedCounts,corelinCountTooverlaps30d6, by = c("cellID", "SampleNum","cellIDSN")) #5860 (4 more than without GFP- cells?)
#pull out unique barcodes and mutate a column to be Bi where i is the row number, then inner_join with logNormalizedCountsSubset
plusnegrepeatBarcodes = plusneglogNormalizedCountsSubsetWBarcodes %>% 
  group_by(BC30StarcodeD6) %>%
  summarize(barCount=n()) #gives number of times each barcode appears - 1024 unique barcodes in 5860 cells (more than without GFP- condition? double check cellID); 1022 unique barcodes in 5860 cells (2/0.3) vs 986 unique barcodes in 5557 cells (2/0.2) vs 504 unique barcodes in 1572 cells vs 515 unique barcodes associated with 2369 unique cells in 3/11 seq run
plusnegreplogNormalizedCountsSubsetWBarcodes = left_join(plusneglogNormalizedCountsSubsetWBarcodes,plusnegrepeatBarcodes, by = "BC30StarcodeD6") %>%
  arrange(desc(barCount)) #5860

uniqueBarcodes <- readRDS(file = paste0(barcodeDirectory,'plusneguniqueBarcodes_30d6.rds')) 
plusneguniqueBarcodes = as_tibble(unique(plusnegreplogNormalizedCountsSubsetWBarcodes$BC30StarcodeD6)) %>%
  transmute(BC30StarcodeD6=value) %>% #can use to load saved barcodes
  inner_join(uniqueBarcodes) #%>%
#mutate(barNum = paste0("B",row_number())) %>%
#transmute(BC30StarcodeD6 = value, barNum)
#save unique barcodes
#saveRDS(plusneguniqueBarcodes, file = paste0(barcodeDirectory,'plusneguniqueBarcodes_30d6.rds'))

#index into umapCoordinates
plusneglinCount_barcodes = inner_join(corelinCountTooverlaps30d6,plusneguniqueBarcodes, by = "BC30StarcodeD6") #5858
plusnegappendumapCoordinates = left_join(plusnegumapCoordinates, plusneglinCount_barcodes, by = c("cellID","SampleNum","cellIDSN")) %>%
  arrange(!is.na(barNum),barNum) #17599

#pull top barcodes
mixedrank = function(x) order(gtools::mixedorder(x))
plusnegbarcodebreakdown = plusneglogNormalizedCountsSubsetWBarcodes %>% 
  group_by(BC30StarcodeD6, SampleNum) %>%
  summarize(barSamCount=n()) %>% #gives number of times each barcode appears per sample number
  spread(SampleNum, barSamCount) %>%
  ungroup() %>%
  inner_join(plusneguniqueBarcodes, by = "BC30StarcodeD6") %>%
  arrange(mixedrank(barNum))

#plotData directory
plotDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Supplement/barcodeB31B_probdistribution/"
KLD <- function(x,y) sum(x * log2(x/y))
JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))

sibBcounts = plusnegreplogNormalizedCountsSubsetWBarcodes %>% #2544 cells
  filter(SampleNum == 'sibB')

sibB_20cellbarcodes = plusnegbarcodebreakdown %>%
  filter(sibB>=20)

sibB_20cellcounts = sibBcounts %>%
  filter(BC30StarcodeD6 %in% sibB_20cellbarcodes$BC30StarcodeD6) %>%
  filter(plusnegcells_clusters!="14")%>%
  select(cellIDSN,plusnegcells_clusters) #950 cells

sibB_20cellbarcodeclusters_new = sibBcounts %>% #select from just one sibling
  inner_join(sibB_20cellbarcodes, by='BC30StarcodeD6') %>% #1578 cells have one of the 20 barcodes
  select(cellIDSN,BC30StarcodeD6,barNum,plusnegcells_clusters,sibB) %>%
  group_by(barNum,BC30StarcodeD6,plusnegcells_clusters, sibB) %>% #group by res=0.5 clusters
  summarise(cellspercluster = n()) %>%
  complete(nesting(barNum,BC30StarcodeD6,sibB),plusnegcells_clusters) %>% #make sure every cluster is listed for each barcode (for graphing purposes only)
  unique() %>% 
  filter(plusnegcells_clusters!="14") %>% #0 cells in cluster 14 -- will throw an error in chi square/Fisher if present (expect value = 0)
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  group_by(plusnegcells_clusters) %>% 
  mutate(allbarperclust = sum(cellspercluster)) %>%
  mutate(propperbarclust = cellspercluster/allbarperclust) %>% #normalize to barcoded cells per cluster (allbarperclust)
  ungroup() %>%
  group_by(barNum,BC30StarcodeD6) %>%
  mutate(prop = propperbarclust/sum(propperbarclust)) %>%
  arrange(mixedrank(barNum)) #14 clusters * 19 barcodes = 266 obs

barperclustB = sibB_20cellbarcodeclusters_new %>%
  ungroup() %>%
  select(plusnegcells_clusters,allbarperclust) %>%
  unique()

barcodeB31B = filter(plusnegbarcodebreakdown,barNum=="B31")
barcode = barcodeB31B$barNum
print(barcode)
bartitle = paste0("barcode ", barcode,"B")
barcodeclusters = sibB_20cellbarcodeclusters_new %>%
  filter(barNum == barcode) #should be 14 rows regardless
allrandomprop = barcodeclusters %>% ungroup() %>% 
  select(plusnegcells_clusters,allbarperclust) %>%
  transmute(plusnegcells_clusters,allbarpercluster=allbarperclust)
allrandomcells = barcodeclusters %>% ungroup() %>%
  select(plusnegcells_clusters)

#randomize
n = 1000 
seeds = 1:n
for (j in seeds){
  set.seed(j) #randomize within cells that are assigned to any of these 30 barcodes
  randomsibBrows = as.integer(sample(rownames(sibB_20cellcounts),barcodeclusters$sibB[1]))
  randomsibB = suppressMessages(sibB_20cellcounts[randomsibBrows,] %>%
                                  select(plusnegcells_clusters) %>%
                                  group_by(plusnegcells_clusters) %>%
                                  summarise(randomcellspercluster = n()) %>%
                                  full_join(barperclustB,by = c('plusnegcells_clusters')) %>%
                                  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
                                  mutate(randomprop = randomcellspercluster/allbarperclust) %>%
                                  mutate(randomprop = randomprop/sum(randomprop)) %>%
                                  arrange(plusnegcells_clusters)) #else zeroes always fall at end
  allrandomprop[,paste0("randomprop",j)] = randomsibB$randomprop
  allrandomcells[,paste0("randomcells",j)] = randomsibB$randomcellspercluster
}
#cell-based
allrandomcells = allrandomcells %>% rowwise() %>%
  mutate(avgrandomcells = mean(c_across(starts_with("randomcells"))))
summaryrandomcells = allrandomcells %>%
  inner_join(barcodeclusters, by='plusnegcells_clusters') %>%
  select(plusnegcells_clusters,cellspercluster,avgrandomcells)

plotsummaryrandomclusters = summaryrandomcells %>%
  transmute(plusnegcells_clusters, observed_cells = cellspercluster, randomavg_cells=avgrandomcells) %>%
  pivot_longer(cols=-plusnegcells_clusters,names_to = c("sample",".value"),names_pattern="(.+)_(.+)")
rawsample.labs <- c(paste0("observed ",barcode,"B"),
                 paste0("average of 1000 random samples"))
names(rawsample.labs) <- c("observed", "randomavg")
rawcellplot = ggplot(plotsummaryrandomclusters, aes(x=plusnegcells_clusters,y=cells,fill=sample)) +
  geom_bar(stat='identity') +
  facet_wrap(~sample,ncol=1, labeller = labeller(sample = rawsample.labs))+
  scale_fill_manual(breaks=c("observed","randomavg"),values = c("#009292","gray80")) +
  ggtitle(paste(bartitle,': ',barcodeclusters$sibB[1],' split B cells\n',sep="")) +
  xlab("Seurat clusters, res=0.5")+
  theme_bw()+
  theme(legend.position = "none")
# save raw cell distribution for B31B in split B
pdf(paste0(plotDirectory,barcode,"_distributionvsavgrandom_rawcells.pdf"), width=4, height=4)
print(rawcellplot)
dev.off()

#proportions
allrandomprop = allrandomprop %>% rowwise() %>%
  mutate(avgrandomprop = mean(c_across(starts_with("randomprop"))))
summaryrandomprop = allrandomprop %>%
  inner_join(barcodeclusters, by='plusnegcells_clusters') %>%
  select(plusnegcells_clusters,prop,avgrandomprop)

#entropy
p_observed = summaryrandomprop %>% ungroup() %>%
  arrange(plusnegcells_clusters)%>%
  summarise(prob_observed = prop+0.000001) #add pseudocount to avoid zero in the numerator and/or denominator of KL
p_random = summaryrandomprop %>% ungroup() %>%
  arrange(plusnegcells_clusters)%>%
  summarise(prob_random = avgrandomprop+0.000001)
  
observedentropy = Entropy(p_observed,base=2) 
randomentropy = Entropy(p_random,base=2) #essentially always 3.708 (i.e. uniform probability x 14 categories)
JS_B <- JSD(p_random,p_observed)
  
plotsummaryrandomprop = summaryrandomprop %>%
  transmute(plusnegcells_clusters, observed_prop = prop, randomavg_prop=avgrandomprop) %>%
  pivot_longer(cols=-plusnegcells_clusters,names_to = c("sample",".value"),names_pattern="(.+)_(.+)")
sample.labs <- c(paste0("observed ",barcode,"B: entropy = ",round(observedentropy,3)),
                 paste0("random sample: entropy = ",round(randomentropy,3)))
names(sample.labs) <- c("observed", "randomavg")
entropyplot_prop = ggplot(plotsummaryrandomprop, aes(x=plusnegcells_clusters,y=prop,fill=sample)) +
  geom_bar(stat='identity') +
  facet_wrap(~sample,ncol=1, labeller = labeller(sample = sample.labs))+
  scale_fill_manual(breaks=c("observed","randomavg"),values = c("#009292","gray80")) +
  ggtitle(paste(bartitle,': ',barcodeclusters$sibB[1],' cells\n','Jensen-Shannon distance = ',round(JS_B,3),sep="")) +
  xlab("Seurat clusters, res=0.5")+
  ylab(paste0("normalized proportion of ",barcode,"B cells"))+
  theme_bw()+
  theme(legend.position = "none")

# save "probability distributions" for B31B in split B
pdf(paste0(plotDirectory,barcode,"_distributionvsavgrandom_NORMno14.pdf"), width=4, height=4)
print(entropyplot_prop)
dev.off()




