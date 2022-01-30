#find actual null number of cells for each cluster
#11/13/21

#pull data from Figure2 (similar to Supp2_barcodeprobabilitydistribution...R)

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
plusneglogNormalizedCountsSubsetWBarcodes = inner_join(plusneglogNormalizedCounts,corelinCountTooverlaps30d6, by = c("cellID", "SampleNum","cellIDSN")) #5860 
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
#mutate(barNum = paste0("B",row_number())) %>% #comment out the above 2 lines and switch with these 2 + saveRDS to save unique barcodes 
#transmute(BC30StarcodeD6 = value, barNum)
#save unique barcodes
#saveRDS(plusneguniqueBarcodes, file = paste0(barcodeDirectory,'plusneguniqueBarcodes_30d6.rds'))

#index into umapCoordinates
plusneglinCount_barcodes = inner_join(corelinCountTooverlaps30d6,plusneguniqueBarcodes, by = "BC30StarcodeD6") #5860
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
plotDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Supplement/"

sibAcounts = plusnegreplogNormalizedCountsSubsetWBarcodes %>% #3316 cells
  filter(SampleNum == 'sibA')

sibA_20cellbarcodes = plusnegbarcodebreakdown %>%
  filter(sibA>=20)

sibA_20cellcounts = sibAcounts %>%
  filter(BC30StarcodeD6 %in% sibA_20cellbarcodes$BC30StarcodeD6) %>% #1578 cells
  filter(plusnegcells_clusters!="14") %>%
  select(cellIDSN,plusnegcells_clusters)#1576 cells: do this bc 14 is so small and sibB has none

sibA_20cellbarcodeclusters_new = sibAcounts %>% #select from just one sibling
  inner_join(sibA_20cellbarcodes, by='BC30StarcodeD6') %>% #1578 cells have one of the 20 barcodes
  select(cellIDSN,BC30StarcodeD6,barNum,plusnegcells_clusters,sibA) %>%
  group_by(barNum,BC30StarcodeD6,plusnegcells_clusters, sibA) %>% #group by res=0.5 clusters
  summarise(cellspercluster = n()) %>%
  complete(nesting(barNum,BC30StarcodeD6,sibA),plusnegcells_clusters) %>% #make sure every cluster is listed for each barcode (for graphing purposes only)
  unique() %>% 
  filter(plusnegcells_clusters!="14") %>% #2 cells in cluster 14 -- remove for sibling B so maybe also remove here for consistency? try it
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  group_by(plusnegcells_clusters) %>% 
  mutate(allbarperclust = sum(cellspercluster)) %>%
  mutate(propperbarclust = cellspercluster/allbarperclust) %>%#normalize to barcoded cells per cluster (allbarperclust)
  ungroup() %>%
  group_by(barNum,BC30StarcodeD6) %>%
  mutate(prop = propperbarclust/sum(propperbarclust)) %>%#normalize within each barcode to add to 1
  arrange(mixedrank(barNum)) #14 clusters * 30 barcodes = 420 obs

barperclustA = sibA_20cellbarcodeclusters_new %>%
  ungroup() %>%
  select(plusnegcells_clusters,allbarperclust) %>%
  unique()

#barcode 28A
barcode28A = filter(plusnegbarcodebreakdown,barNum=="B28")

barcode = barcode28A$barNum
print(barcode)
bartitle = paste0("barcode ", barcode,"A")
barcodeclusters = sibA_20cellbarcodeclusters_new %>%
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
  randomsibArows = as.integer(sample(rownames(sibA_20cellcounts),barcodeclusters$sibA[1]))
  randomsibA = suppressMessages(sibA_20cellcounts[randomsibArows,] %>%
                                  select(plusnegcells_clusters) %>%
                                  group_by(plusnegcells_clusters) %>%
                                  summarise(randomcellspercluster = n()) %>%
                                  full_join(barperclustA,by = c('plusnegcells_clusters')) %>%
                                  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
                                  mutate(randomprop = randomcellspercluster/allbarperclust) %>%
                                  mutate(randomprop = randomprop/sum(randomprop)) %>%
                                  arrange(plusnegcells_clusters)) #else zeroes always fall at end
  allrandomprop[,paste0("randomprop",j)] = randomsibA$randomprop
  allrandomcells[,paste0("randomcells",j)] = randomsibA$randomcellspercluster
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
rawsample.labs <- c(paste0("observed ",barcode,"A"),
                    paste0("average of 1000 random samples"))
names(rawsample.labs) <- c("observed", "randomavg")
rawcellplot = ggplot(plotsummaryrandomclusters, aes(x=plusnegcells_clusters,y=cells,fill=sample)) +
  geom_bar(stat='identity') +
  facet_wrap(~sample,ncol=1, labeller = labeller(sample = rawsample.labs))+
  scale_fill_manual(breaks=c("observed","randomavg"),values = c("#009292","gray80")) +
  ggtitle(paste(bartitle,': ',barcodeclusters$sibA[1],' split A cells\n',sep="")) +
  xlab("Seurat clusters, res=0.5")+
  theme_bw()+
  theme(legend.position = "none")

#what is the chance of getting ≥2 cells in cluster 3?
#out of 1000 samples, pull cluster 3 info
randomcluster3cells = as.numeric(allrandomcells %>%
                                   filter(plusnegcells_clusters == "3") %>%
                                   select(starts_with("randomcells")))
#there are 2 cells observed in cluster 3 for barcode 28A
observedcluster3cells = as.numeric(summaryrandomcells %>%
                                     filter(plusnegcells_clusters == "3") %>%
                                     select(cellspercluster))
#distribution is not normal; most trials have 0 in cluster 3 --> right-skewed
#calculate area under the curve
p = sum(randomcluster3cells>=observedcluster3cells)/n

B28Acluster3significance = ggplot() +
  geom_histogram(aes(x=randomcluster3cells),binwidth=1,fill='gray80') +
  geom_vline(xintercept=observedcluster3cells,color='#009292') +
  ggtitle(paste0('barcode 28A cluster 3 \nP(≥',observedcluster3cells,' cells) = ',p)) +
  xlab("number of cells in cluster 3") +
  #coord_cartesian(xlim = c(0, 1))+ #without coord_cartesian throws an error about rows being removed
  theme_bw()
#save plot showing observed JSD vs 1000 random samples
cairo_pdf(paste0(plotDirectory, barcode,"Acluster3distribution.pdf"), width=4, height=4) #allows use of ≥
print(B28Acluster3significance)
dev.off()

#barcode 35A
barcode35A = filter(plusnegbarcodebreakdown,barNum=="B35")

barcode = barcode35A$barNum
print(barcode)
bartitle = paste0("barcode ", barcode,"A")
barcodeclusters = sibA_20cellbarcodeclusters_new %>%
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
  randomsibArows = as.integer(sample(rownames(sibA_20cellcounts),barcodeclusters$sibA[1]))
  randomsibA = suppressMessages(sibA_20cellcounts[randomsibArows,] %>%
                                  select(plusnegcells_clusters) %>%
                                  group_by(plusnegcells_clusters) %>%
                                  summarise(randomcellspercluster = n()) %>%
                                  full_join(barperclustA,by = c('plusnegcells_clusters')) %>%
                                  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
                                  mutate(randomprop = randomcellspercluster/allbarperclust) %>%
                                  mutate(randomprop = randomprop/sum(randomprop)) %>%
                                  arrange(plusnegcells_clusters)) #else zeroes always fall at end
  allrandomprop[,paste0("randomprop",j)] = randomsibA$randomprop
  allrandomcells[,paste0("randomcells",j)] = randomsibA$randomcellspercluster
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
rawsample.labs <- c(paste0("observed ",barcode,"A"),
                    paste0("average of 1000 random samples"))
names(rawsample.labs) <- c("observed", "randomavg")
rawcellplot = ggplot(plotsummaryrandomclusters, aes(x=plusnegcells_clusters,y=cells,fill=sample)) +
  geom_bar(stat='identity') +
  facet_wrap(~sample,ncol=1, labeller = labeller(sample = rawsample.labs))+
  scale_fill_manual(breaks=c("observed","randomavg"),values = c("#009292","gray80")) +
  ggtitle(paste(bartitle,': ',barcodeclusters$sibA[1],' split A cells\n',sep="")) +
  xlab("Seurat clusters, res=0.5")+
  theme_bw()+
  theme(legend.position = "none")

#what is the chance of getting ≥2 cells in cluster 3?
#out of 1000 samples, pull cluster 3 info
randomcluster3cells = as.numeric(allrandomcells %>%
                                   filter(plusnegcells_clusters == "3") %>%
                                   select(starts_with("randomcells")))
#there are 2 cells observed in cluster 3 for barcode 35A
observedcluster3cells = as.numeric(summaryrandomcells %>%
                                     filter(plusnegcells_clusters == "3") %>%
                                     select(cellspercluster))
#distribution is not normal; most trials have 0 in cluster 3 --> right-skewed
#calculate area under the curve
p = sum(randomcluster3cells>=observedcluster3cells)/n

B35Acluster3significance = ggplot() +
  geom_histogram(aes(x=randomcluster3cells),binwidth=1,fill='gray80') +
  geom_vline(xintercept=observedcluster3cells,color='#009292') +
  ggtitle(paste0('barcode 35A cluster 3 \nP(≥',observedcluster3cells,' cells) = ',p)) +
  xlab("number of cells in cluster 3") +
  #coord_cartesian(xlim = c(0, 1))+ #without coord_cartesian throws an error about rows being removed
  theme_bw()
#save plot showing observed JSD vs 1000 random samples
cairo_pdf(paste0(plotDirectory, barcode,"Acluster3distribution.pdf"), width=4, height=4) #allows use of ≥
print(B35Acluster3significance)
dev.off()