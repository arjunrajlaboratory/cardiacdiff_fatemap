#Figure 2: Constraint analysis
#July2021

library(tidyverse)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(grid)
library(DescTools)
###############Function for ggplot by Lee Richman
create_lpr_theme <- function(){
  lpr_theme <- ggplot2::theme_bw() + ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 32),
                   plot.title = ggplot2::element_text(face = "bold",
                                                      hjust = 0.5, size = ggplot2::rel(0.5)), axis.text = ggplot2::element_text(size = ggplot2::rel(0.8),
                                                                                                                                color = "black", face = "bold"), axis.title = ggplot2::element_text(size = ggplot2::rel(0.6),
                                                                                                                                                                                                    face = "bold"), legend.title = ggplot2::element_blank(),
                   legend.position = "bottom", axis.line = element_line(size = 2),
                   axis.ticks = element_line(size = 2), strip.background = element_rect(size = 2),
                   strip.text = element_text(size = ggplot2::rel(0.7),
                                             face = "bold"))
  return(lpr_theme)
}
##################

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


#Figure 2: 
#plotData directory
plotDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Figure2/"
#Jensen-Shannon divergence: use https://enterotype.embl.de/enterotypes.html
KLD <- function(x,y) sum(x * log2(x/y))
JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))

#Figure 2A - only keep cells with a barcode for UMAP schematic
barcodedcellsUMAP <- ggplot() +
  geom_point(data = plusnegappendumapCoordinates[which(plusnegappendumapCoordinates$barNum !="NA"),], aes(x = UMAP_1, y = UMAP_2, color = "gray95"), size = 1) +
  scale_color_manual(values = c("gray80"), labels = c("barcoded")) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.3)),
        legend.text = element_text(size = rel(0.3), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())

pdf(paste(plotDirectory,"UMAPschematics/barcodeonlyUMAP.pdf", sep=""), width=4.55, height=4.2)
print(barcodedcellsUMAP)
dev.off()


#Figure 2B: run JSD analysis on 30 barcodes in split A and 19 barcodes in split B that label ≥20 cells
###############split A analysis
sibAcounts = plusnegreplogNormalizedCountsSubsetWBarcodes %>% #3316 cells
  filter(SampleNum == 'sibA')

#30 barcodes with more than 20 sibling A cells (without considering number of sibling B cells)
sibA_20cellbarcodes = plusnegbarcodebreakdown %>%
  filter(sibA>=20)
#"B1"  "B2"  "B3"  "B4"  "B5"  "B6"  "B8"  "B9"  "B10" "B11" 
#"B12" "B13" "B14" "B17" "B18" "B19" "B22" "B25" "B28" "B29" 
#"B30" "B34" "B35" "B36" "B39" "B44" "B45" "B46" "B49" "B51"

sibA_20cellcounts = sibAcounts %>%
  filter(BC30StarcodeD6 %in% sibA_20cellbarcodes$BC30StarcodeD6) %>% #1578 cells
  filter(plusnegcells_clusters!="14") %>%
  select(cellIDSN,plusnegcells_clusters)#1576 cells: do this bc 14 is so small and sibB has none

sibA_20cellumap = plusnegappendumapCoordinates[which((plusnegappendumapCoordinates$barNum %in% sibA_20cellbarcodes$barNum)&(plusnegappendumapCoordinates$SampleNum=="sibA")),]

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

#save RDS file that is used in Supp2_constraintheatmap.R script
#MAKE SURE TO RUN THIS BEFORE RUNNING THAT SCRIPT
saveRDS(sibA_20cellbarcodeclusters_new, file = paste0(plotDirectory,'sibA_20cellbarcodeclusters_new.rds'))

barperclustA = sibA_20cellbarcodeclusters_new %>%
  ungroup() %>%
  select(plusnegcells_clusters,allbarperclust) %>%
  unique()

n = 1000 
seeds = 1:n
JStest = rep(0,n)
JS_A=rep(0,length(sibA_20cellbarcodes$barNum))
meanJS_A=rep(0,length(sibA_20cellbarcodes$barNum))
sdJS_A=rep(0,length(sibA_20cellbarcodes$barNum))
setwd(paste0(plotDirectory,"sibAconstraint_res0.5/"))
for (i in 1:length(sibA_20cellbarcodes$barNum)){
  barcode = sibA_20cellbarcodes$barNum[i]
  print(barcode)
  bartitle = paste0("barcode ", barcode,"A")
  barcodeclusters = sibA_20cellbarcodeclusters_new %>%
    filter(barNum == barcode) #should be 14 rows regardless
  allrandomprop = barperclustA
  #randomize
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
  }
  #for entropy analysis
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
  JS_A[i] <- JSD(p_random,p_observed)
  
  plotsummaryrandomprop = summaryrandomprop %>%
    transmute(plusnegcells_clusters, observed_prop = prop, randomavg_prop=avgrandomprop) %>%
    pivot_longer(cols=-plusnegcells_clusters,names_to = c("sample",".value"),names_pattern="(.+)_(.+)")
  sample.labs <- c(paste0("observed ",barcode,"A: entropy = ",round(observedentropy,3)),
                   paste0("random sample: entropy = ",round(randomentropy,3)))
  names(sample.labs) <- c("observed", "randomavg")
  entropyplot_prop = ggplot(plotsummaryrandomprop, aes(x=plusnegcells_clusters,y=prop,fill=sample)) +
    geom_bar(stat='identity') +
    facet_wrap(~sample,ncol=1, labeller = labeller(sample = sample.labs))+
    scale_fill_manual(breaks=c("observed","randomavg"),values = c("#009292","gray80")) +
    ggtitle(paste(bartitle,': ',barcodeclusters$sibA[1],' cells\n','Jensen-Shannon distance = ',round(JS_A[i],3),sep="")) +
    xlab("Seurat clusters, res=0.5")+
    ylab(paste0("normalized proportion of ",barcode,"A cells"))+
    theme_bw()+
    theme(legend.position = "none")
  # save "probability distributions" for each barcode in split A
  pdf(paste(barcode,"distributionvsavgrandom_NORMno14.pdf", sep="_"), width=4, height=4)
  print(entropyplot_prop)
  dev.off()
  
  #test JSD significance
  for (k in seeds+n){ #seeds+n will give 1001-2000
    set.seed(k)
    JSDsibArows = as.integer(sample(rownames(sibA_20cellcounts),barcodeclusters$sibA[1]))
    JSDsibA = sibA_20cellcounts[JSDsibArows,] %>%
      select(plusnegcells_clusters) %>%
      full_join(barperclustA,by = 'plusnegcells_clusters') %>%
      group_by(plusnegcells_clusters,allbarperclust) %>%
      summarise(cellspercluster = n(),.groups='drop_last') %>% #drop_last > keep bc result is the same
      mutate(propperbarclust = cellspercluster/allbarperclust) %>%
      ungroup() %>%
      mutate(prop = propperbarclust/sum(propperbarclust))
    p_test = JSDsibA %>% ungroup() %>%
      arrange(plusnegcells_clusters)%>%
      summarise(prob_test = prop+0.000001)
    JStest[k-n] <- JSD(p_random,p_test)
  }
  meanJS_A[i] = mean(JStest)
  sdJS_A[i] = sd(JStest)
  normp = pnorm(JS_A[i],mean=meanJS_A[i],sd=sdJS_A[i],lower.tail=FALSE) #have to edit to get upper tail?
  if(round(normp,3)==0){
    normp3 <- "< 0.001"
  }else{
    normp3 <- paste("=",round(normp,3))
  }
  JSplot = ggplot() +
    geom_histogram(aes(x=JStest),binwidth=0.01,fill='gray80') +
    geom_vline(xintercept=JS_A[i],color='#009292') +
    ggtitle(paste0(bartitle,': ',barcodeclusters$sibA[1],' cells','\np ',normp3)) +
    xlab("Jensen-Shannon distance") +
    coord_cartesian(xlim = c(0, 1))+ #without coord_cartesian throws an error about rows being removed
    theme_bw()
  #save plot showing observed JSD vs 1000 random samples
  pdf(paste(barcode,"JSDdistribution.pdf", sep="_"), width=4, height=4)
  print(JSplot)
  dev.off()
  
  #walk through different barcode in split A only
  umap = ggplot() +
    geom_point(data = sibA_20cellumap, aes(x = UMAP_1, y = UMAP_2, color = "other")) +
    #geom_point(data = sibA_20cellumap[which(sibA_20cellumap$plusnegcells_clusters%in%top3clusters$plusnegcells_clusters),],aes(x = UMAP_1, y = UMAP_2, color = "top3clust"),size=1)+
    geom_point(data = sibA_20cellumap[which((sibA_20cellumap$barNum==barcode)),], aes(x = UMAP_1, y = UMAP_2, color = "barcode"), size = 3) +
    scale_color_manual(breaks=c("barcode","other"),values = c("#009292","gray80"), labels = c(paste0(barcode),"other")) +
    create_lpr_theme() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = rel(0.3)),
          legend.text = element_text(size = rel(0.3), angle = 30),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    labs(color =  'Split A barcodes\nwith at least 20 cells') +
    ggtitle(paste0(bartitle))
  # save UMAP showing barcoded split A cells in teal
  pdf(paste(barcode,"umap_alone.pdf", sep="_"), width=6, height=6)
  print(umap)
  dev.off()
}
#save JSD values to create table
JStable_A = tibble(JS_A,meanJS_A,sdJS_A,barNum=sibA_20cellbarcodes$barNum,nCells=sibA_20cellbarcodes$sibA)
write_csv(JStable_A, file = paste0(plotDirectory,"JStableinfo_sibA.csv"))

###############split B alone
sibBcounts = plusnegreplogNormalizedCountsSubsetWBarcodes %>% #2544 cells
  filter(SampleNum == 'sibB')

#19 barcodes with more than 20 sibling B cells (without considering number of sibling A cells)
sibB_20cellbarcodes = plusnegbarcodebreakdown %>%
  filter(sibB>=20)
#"B1"  "B2"  "B4"  "B7"  "B15" "B16" "B20" "B24" "B26" "B27" 
#"B31" "B32" "B33" "B38" "B40" "B41" "B42" "B43" "B48" 

sibB_20cellcounts = sibBcounts %>%
  filter(BC30StarcodeD6 %in% sibB_20cellbarcodes$BC30StarcodeD6) %>%
  filter(plusnegcells_clusters!="14")%>%
  select(cellIDSN,plusnegcells_clusters) #950 cells

sibB_20cellumap = plusnegappendumapCoordinates[which((plusnegappendumapCoordinates$barNum %in% sibB_20cellbarcodes$barNum)&(plusnegappendumapCoordinates$SampleNum=="sibB")),]

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

#save RDS file that is used in Supp2_constraintheatmap.R script
#MAKE SURE TO RUN THIS BEFORE RUNNING THAT SCRIPT
saveRDS(sibB_20cellbarcodeclusters_new, file = paste0(plotDirectory,'sibB_20cellbarcodeclusters_new.rds'))

barperclustB = sibB_20cellbarcodeclusters_new %>%
  ungroup() %>%
  select(plusnegcells_clusters,allbarperclust) %>%
  unique()

JS_B=rep(0,length(sibB_20cellbarcodes$barNum))
meanJS_B=rep(0,length(sibB_20cellbarcodes$barNum))
sdJS_B=rep(0,length(sibB_20cellbarcodes$barNum))
setwd(paste0(plotDirectory,"sibBconstraint_res0.5/"))
for (i in 1:length(sibB_20cellbarcodes$barNum)){
  barcode = sibB_20cellbarcodes$barNum[i]
  print(barcode)
  bartitle = paste0("barcode ", barcode,"B")
  barcodeclusters = sibB_20cellbarcodeclusters_new %>%
    filter(barNum == barcode) #should be 14 rows regardless
  allrandomprop = barcodeclusters %>% ungroup() %>% select(plusnegcells_clusters,allbarperclust) %>%
    transmute(plusnegcells_clusters,allbarpercluster=allbarperclust)
  #randomize
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
  }
  #for entropy analysis
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
  JS_B[i] <- JSD(p_random,p_observed)
  
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
    ggtitle(paste(bartitle,': ',barcodeclusters$sibB[1],' cells\n','Jensen-Shannon distance = ',round(JS_B[i],3),sep="")) +
    xlab("Seurat clusters, res=0.5")+
    ylab(paste0("normalized proportion of ",barcode,"B cells"))+
    #scale_y_continuous(breaks=seq(0,0.2,0.1),limits=c(0,0.23))+ #edited B16B to be better for AI
    theme_bw()+
    theme(legend.position = "none")
  # save "probability distributions" for each barcode in split B
  pdf(paste(barcode,"distributionvsavgrandom_NORMno14.pdf", sep="_"), width=4, height=4)
  print(entropyplot_prop)
  dev.off()
  
  #test JS test significance
  for (k in seeds+n){ #seeds+n will give 1001-2000
    set.seed(k)
    JSDsibBrows = as.integer(sample(rownames(sibB_20cellcounts),barcodeclusters$sibB[1]))
    JSDsibB = sibB_20cellcounts[JSDsibBrows,] %>%
      select(plusnegcells_clusters) %>%
      full_join(barperclustB,by = 'plusnegcells_clusters') %>%
      group_by(plusnegcells_clusters,allbarperclust) %>%
      summarise(cellspercluster = n(),.groups='drop_last') %>% #drop_last > keep bc result is the same
      mutate(propperbarclust = cellspercluster/allbarperclust) %>%
      ungroup() %>%
      mutate(prop = propperbarclust/sum(propperbarclust))
    p_test = JSDsibB %>% ungroup() %>%
      arrange(plusnegcells_clusters)%>%
      summarise(prob_test = prop+0.000001)
    JStest[k-n] <- JSD(p_random,p_test)
  }
  meanJS_B[i] = mean(JStest)
  sdJS_B[i] = sd(JStest)
  normp = pnorm(JS_B[i],mean=meanJS_B[i],sd=sdJS_B[i],lower.tail=FALSE) #have to edit to get upper tail?
  if(round(normp,3)==0){
    normp3 <- "< 0.001"
  }else{
    normp3 <- paste("=",round(normp,3))
  }
  JSplot = ggplot() +
    geom_histogram(aes(x=JStest),binwidth=0.01,fill='gray80') +
    geom_vline(xintercept=JS_B[i],color='#009292') +
    ggtitle(paste0(bartitle,': ',barcodeclusters$sibB[1],' cells','\np ',normp3)) +
    xlab("Jensen-Shannon distance") +
    coord_cartesian(xlim = c(0, 1))+ #without coord_cartesian throws an error about rows being removed
    theme_bw()
  #save plot showing observed JSD vs 1000 random samples
  pdf(paste(barcode,"JSDdistribution.pdf", sep="_"), width=4, height=4)
  print(JSplot)
  dev.off()
  
  #walk through different barcode in split B only
  umap = ggplot() +
    geom_point(data = sibB_20cellumap, aes(x = UMAP_1, y = UMAP_2, color = "other")) +
    #geom_point(data = sibB_20cellumap[which(sibB_20cellumap$plusnegcells_clusters%in%top3clusters$plusnegcells_clusters),],aes(x = UMAP_1, y = UMAP_2, color = "top3clust"),size=1)+
    geom_point(data = sibB_20cellumap[which((sibB_20cellumap$barNum==barcode)),], aes(x = UMAP_1, y = UMAP_2, color = "barcode"), size = 3) +
    scale_color_manual(breaks=c("barcode","other"),values = c("#009292","gray80"), labels = c(paste0(barcode),"other")) +
    create_lpr_theme() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = rel(0.3)),
          legend.text = element_text(size = rel(0.3), angle = 30),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    labs(color =  'Split B barcodes\nwith at least 20 cells') +
    ggtitle(paste0(bartitle))
  # save UMAP showing barcoded split B cells in teal
  pdf(paste(barcode,"umap_alone.pdf", sep="_"), width=6, height=6)
  print(umap)
  dev.off()
}
#save JSD values to create table
JStable_B = tibble(JS_B,meanJS_B,sdJS_B,barNum=sibB_20cellbarcodes$barNum,nCells=sibB_20cellbarcodes$sibB)
write_csv(JStable_B, file = paste0(plotDirectory,"JStableinfo_sibB.csv"))

