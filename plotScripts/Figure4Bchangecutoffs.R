#Rerun figure 4B with different cell number cutoffs
#10/31/21

library(tidyverse)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
#load Seurat and barcode data as before
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
plotDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Supplement/memory/"
#Jensen-Shannon divergence: use https://enterotype.embl.de/enterotypes.html
KLD <- function(x,y) sum(x * log2(x/y))
JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))

###############pull all barcodes with ≥N cells per split
plusnegbarcodebreakdown_noNA = plusnegbarcodebreakdown %>% na.omit() #81 of 1022 (2/0.3)
N = c(7, 10, 15, 20)
for (n in 1:length(N)){
  plusnegbarcodebreakdown_noNA_more = plusnegbarcodebreakdown_noNA %>% 
    filter(sibA>=N[n]) %>% 
    filter(sibB>=N[n]) #18 barcodes for 5

  cellsperclustperSampleNum <- plusneglogNormalizedCountsSubsetWBarcodes %>% 
    group_by(BC30StarcodeD6, SampleNum, plusnegcells_clusters) %>% #group by res=0.5 clusters
    summarize(barSamclustCount=n()) %>% #gives number of cells per barcode per sibling per cluster
    inner_join(plusnegbarcodebreakdown_noNA_more, by = "BC30StarcodeD6") %>%
    select(-sibA,-sibB) %>% 
    complete(barNum,plusnegcells_clusters,SampleNum) %>%
    mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
    filter(plusnegcells_clusters!="14") %>% #0 cells in cluster 14 - will throw an error if not removed
    group_by(SampleNum,plusnegcells_clusters) %>%
    mutate(allbarperclust = sum(barSamclustCount)) %>%
    mutate(propperbarclust = barSamclustCount/allbarperclust) %>%
    ungroup() %>%
    group_by(barNum,BC30StarcodeD6,SampleNum) %>%
    mutate(prop = propperbarclust/sum(propperbarclust)) %>%
    arrange(mixedrank(barNum)) 
  
  #heatmap of pairwise JSD
  propsibA = cellsperclustperSampleNum %>%
    filter(SampleNum == 'sibA') %>%
    ungroup() %>%
    mutate(prob = prop+0.000001) %>% #add to avoid NaNs due to zeroes
    select(barNum,plusnegcells_clusters,prob) %>%
    pivot_wider(id_cols = plusnegcells_clusters,names_prefix = "probA",names_from = barNum,values_from = prob)
  
  propsibB = cellsperclustperSampleNum %>%
    filter(SampleNum == 'sibB') %>%
    ungroup() %>%
    mutate(prob = prop+0.000001) %>% #add to avoid NaNs due to zeroes
    select(barNum,plusnegcells_clusters,prob) %>%
    pivot_wider(id_cols = plusnegcells_clusters,names_prefix = "probB",names_from = barNum,values_from = prob)
  
  JSDMatrix = matrix(, nrow = nrow(plusnegbarcodebreakdown_noNA_more), ncol = nrow(plusnegbarcodebreakdown_noNA_more))
  #generates empty matrix that is 18x18 (for the 18 barcodes that are associated with at least 5 cells in both arms)
  for (i in c(2:length(propsibA))) {
    for (j in c(2:length(propsibB))) {
      JSDMatrix[i-1,j-1] = JSD(propsibA[,i],propsibB[,j]) #start at 2 so have to index at i-1, j-1
    }
  }
  
  usenames = plusnegbarcodebreakdown_noNA_more$barNum
  JSDMatrixplot = JSDMatrix
  rownames(JSDMatrixplot) = sub("$", "A", usenames)#colnames(propsibA)[2:19]
  colnames(JSDMatrixplot) = sub("$", "B", usenames)#colnames(propsibB)[2:19]
  melted_JSDMatrix <- melt(JSDMatrixplot, na.rm = TRUE) #list of JSD values
  melted_JSDMatrix$value <- round(melted_JSDMatrix$value,1) #round to 1 decimal point
  #reverse order per https://stackoverflow.com/questions/40725146/how-do-i-reverse-the-y-axis-in-ggplots-geom-tile
  
  #pull for outlined rectangles
  diag_JSDMatrix <- melted_JSDMatrix[which(melted_JSDMatrix$Var2 == sub("A", "B", melted_JSDMatrix$Var1)),]
  plots_JSDMatrix <- na.omit(diag_JSDMatrix[which(melted_JSDMatrix$Var1 %in% c("B1A","B2A","B14A","B23A")),])
  
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
    geom_rect(data=diag_JSDMatrix, size=1, fill=NA, colour="black", #used to use #FFFF6D",
              aes(xmin=as.integer(Var2) - 0.5, xmax=as.integer(Var2) + 0.5, ymin=as.integer(forcats::fct_rev(Var1)) - 0.5, ymax=as.integer(forcats::fct_rev(Var1)) + 0.5)) +
    geom_rect(data=plots_JSDMatrix, size=1, fill=NA, colour="#ff068f",
              aes(xmin=as.integer(Var2) - 0.5, xmax=as.integer(Var2) + 0.5, ymin=as.integer(forcats::fct_rev(Var1)) - 0.5, ymax=as.integer(forcats::fct_rev(Var1)) + 0.5)) +
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
      legend.position = "bottom",
      legend.direction = "horizontal",
      text=element_text(family="Helvetica"))+
    guides(fill = guide_colorbar(barwidth = 9, barheight = 1,
                                 title.position = "top", title.hjust = 0.5)) +
    scale_y_discrete(position = "right")
  
  #save JSD pairwise comparison heatmap
  pdf(paste0(plotDirectory,N[n],"cells/JSDheatmap_1dec.pdf"), width=7, height=7)
  print(plot_JSD)
  dev.off()
  
  
  #barplots and umap
  setwd(paste0(plotDirectory))
  for (i in c(1:nrow(plusnegbarcodebreakdown_noNA_more))){
    barcode = plusnegbarcodebreakdown_noNA_more$barNum[i]; print(barcode)
    barclusters = cellsperclustperSampleNum %>%
      filter(barNum == barcode) #%>% #should be 14 rows regardless
    bartitle = paste("barcode", barcode)
    
    #create probability distributions
    p_sibA = barclusters %>%
      filter(SampleNum == 'sibA') %>%
      ungroup() %>%
      arrange(plusnegcells_clusters)%>%
      summarise(prob_A = prop+0.000001)
    p_sibB = barclusters %>%
      filter(SampleNum == 'sibB') %>%
      ungroup() %>%
      arrange(plusnegcells_clusters)%>%
      summarise(prob_B = prop+0.000001)
    
    JS <- JSD(p_sibB,p_sibA)
    
    samplenum.labs <- c(paste0("barcode ",barcode,"A",": ",plusnegbarcodebreakdown_noNA_more[i,]$sibA," cells"),
                        paste0("barcode ",barcode,"B",": ",plusnegbarcodebreakdown_noNA_more[i,]$sibB," cells"))
    names(samplenum.labs) <- c("sibA", "sibB")
    plotnormcellsclustsample <- ggplot(barclusters,aes(x=plusnegcells_clusters,y=prop,group=SampleNum,fill=SampleNum)) +
      geom_bar(stat='identity')+
      facet_wrap(~SampleNum,ncol=1,labeller = labeller(SampleNum = samplenum.labs))+
      ylab(paste0("normalized proportion of ",barcode," cells"))+xlab("Seurat clusters, res=0.5")+
      scale_fill_manual(breaks=c("sibA","sibB"),values = c("#FFB6DB","#ff068f")) +
      ggtitle(paste(barcode,"A vs ", barcode,"B\nJensen-Shannon distance = ",round(JS,3),sep=""))+
      theme_bw()+
      theme(legend.position = "none")
    
    #save "probability distribution" bargraphs for each barcode
    pdf(paste0(plotDirectory,N[n],"cells/",barcode,"normcellsperclusterpersample.pdf"), width=4, height=4)
    print(plotnormcellsclustsample)
    dev.off()
    
    umap = ggplot() +
      geom_point(data = plusnegappendumapCoordinates[which((plusnegappendumapCoordinates$barNum %in% plusnegbarcodebreakdown_noNA_more$barNum)),], aes(x = UMAP_1, y = UMAP_2, color = "other"), size = 1)+
      geom_point(data = plusnegappendumapCoordinates[which((plusnegappendumapCoordinates$barNum==barcode) & (plusnegappendumapCoordinates$SampleNum=="sibA")),], aes(x = UMAP_1, y = UMAP_2, color = "sibA"), size = 3) +
      geom_point(data = plusnegappendumapCoordinates[which((plusnegappendumapCoordinates$barNum==barcode) & (plusnegappendumapCoordinates$SampleNum=="sibB")),], aes(x = UMAP_1, y = UMAP_2, color = "sibB"), size = 3) +
      scale_color_manual(breaks = c("sibA","sibB","other"), values = c("sibA"="#FFB6DB","sibB"="#ff068f","other"="gray80"), labels = c(paste("barcode ",barcode,"A",sep=""),paste("barcode ",barcode,"B",sep=""),"other cells")) +
      create_lpr_theme() +
      theme(legend.position = "bottom",
            legend.title = element_text(size = rel(0.3)),
            legend.text = element_text(size = rel(0.3), angle = 20),
            axis.text = element_blank(),
            axis.ticks = element_blank()) +
      labs(color =  "Sample origin of barcoded cells") +
      ggtitle(paste0(bartitle))
    
    #save UMAP marking barcoded cells by which split they were found in
    pdf(paste0(plotDirectory,N[n],"cells/",barcode,"_umap.pdf"), width=6, height=6)
    print(umap)
    dev.off()
  }
  
}

