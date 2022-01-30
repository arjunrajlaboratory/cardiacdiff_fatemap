#Figure 3: survival, top100 and pullouts analysis
#July2021

library(tidyverse)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(grid)
library(VennDiagram)
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

#plotData directory
plotDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Figure3/"

#Figure 3B: survival analysis
#pull barcode info
barcodeDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/extractedData/gDNAbarcodeobjects/"
sampleTable_g1_MOI0.5 <- readRDS(paste0(barcodeDirectory,'sampleTable_g1_MOI0.5.rds'))

####split barcodes into iPS and differentiated objects
sampleTable_g1_MOI0.5_30d8_diff = sampleTable_g1_MOI0.5 %>%
  filter(SampleNum %in% c("diff_05A_g1","diff_05B_g1")) %>%
  dplyr::select(countsperUMI,BC30StarcodeD8,SampleNum) %>% #293974 UMIs for A, 269517 for B, using 30d8
  group_by(BC30StarcodeD8, SampleNum) %>%
  summarise(UMI = n(), counts = sum(countsperUMI)) #collapse and count how many distinct UMIs per barcode
#142854 barcodes bc not collapsing between samples

sampleTable_g1_MOI0.5_30d8_iPS = sampleTable_g1_MOI0.5 %>%
  filter(SampleNum %in% c("iPS_05_g1")) %>%
  dplyr::select(countsperUMI,BC30StarcodeD8,SampleNum) %>% #4478530 UMIs using 30d8
  group_by(BC30StarcodeD8, SampleNum) %>%
  summarise(UMI = n(), counts = sum(countsperUMI)) #collapse and count how many distinct UMIs per barcode
#111649 barcodes - approx 2x expected

####determine actual observed overlap
splitAbar = sampleTable_g1_MOI0.5_30d8_diff %>%
  filter(SampleNum == "diff_05A_g1") %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 72782 barcodes -- a bit more than expected w/o selection
splitBbar = sampleTable_g1_MOI0.5_30d8_diff %>%
  filter(SampleNum == "diff_05B_g1") %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 70072 barcodes
split_initialbar = sampleTable_g1_MOI0.5_30d8_iPS %>%
  filter(SampleNum == "iPS_05_g1")%>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 111649 barcodes -- about 2x expected
venn.diagram(
  x = list(unlist(splitAbar), unlist(splitBbar)),
  category.names = c("splitA" , "splitB"),
  filename = paste0(plotDirectory,'VD_gDNA/g1MOI0.5/nocutoff/g1MOI0.5_VDdiffs.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-15, 15),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
) #more overlap than for MOI 0.1??
venn.diagram(
  x = list(unlist(split_initialbar), unlist(splitAbar)),
  category.names = c("iPS" , "splitA"),
  filename =  paste0(plotDirectory,'VD_gDNA/g1MOI0.5/nocutoff/g1MOI0.5_VDsplitAiPS.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)
venn.diagram(
  x = list(unlist(split_initialbar), unlist(splitBbar)),
  category.names = c("iPS" , "splitB"),
  filename =  paste0(plotDirectory,'VD_gDNA/g1MOI0.5/nocutoff/g1MOI0.5_VDsplitBiPS.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)

####run simulations
alliPS = sampleTable_g1_MOI0.5_30d8_iPS %>%
  filter(SampleNum %in% c("iPS_05_g1")) %>%
  dplyr::select(BC30StarcodeD8) 

#set up null survivability of iPS -- can skip this section after running once
#NEED TO RUN THIS SCRIPT BEFORE RUNNING Supp3_survivaloverlapsignificance.R
survive <- function(datarow){ #reads in a row from sampleTable_g2_MOI1.0_30d8_iPS, designates binomial split by number of UMIs
  barcode_binom = rbinom(datarow[3],1,0.5)
}
#for each barcode, split UMIs across 2 splits
splitA = vector()
splitB = vector()
set.seed(2059)
binomsplit = apply(sampleTable_g1_MOI0.5_30d8_iPS,1,survive)
for (i in 1:nrow(sampleTable_g1_MOI0.5_30d8_iPS)){
  barcode = sampleTable_g1_MOI0.5_30d8_iPS$BC30StarcodeD8[i]
  splitA = c(splitA, rep(barcode,sum(binomsplit[[i]]==0)))
  splitB = c(splitB, rep(barcode,sum(binomsplit[[i]]==1)))
} #took ~2-3 hours
saveRDS(splitA,file=paste0(plotDirectory,'VD_gDNA/g1MOI0.5/simulation/splitA.rds'))
saveRDS(splitB,file=paste0(plotDirectory,'VD_gDNA/g1MOI0.5/simulation/splitB.rds'))

#skip to here to read in RDS files
splitA = readRDS(file=paste0(plotDirectory,'VD_gDNA/g1MOI0.5/simulation/splitA.rds'))
splitB = readRDS(file=paste0(plotDirectory,'VD_gDNA/g1MOI0.5/simulation/splitB.rds'))

survivability = seq(0.06, 0.15, 0.01) #toggle survivability to factor in loss of some UMIs
set.seed(2059)
for (s in 1:length(survivability)){
  survivalA <- sample(length(splitA),length(splitA)*survivability[s])
  splitAsurvive <- splitA[survivalA]
  survivalB <- sample(length(splitB),length(splitB)*survivability[s])
  splitBsurvive <- splitB[survivalB]
  venn.diagram(
    x = list(unlist(alliPS), splitAsurvive),
    category.names = c("iPS" , "splitAsim"),
    filename = paste0(plotDirectory,'VD_gDNA/g1MOI0.5/simulation/splitAiPS_survivability',survivability[s],'.pdf'),
    output=TRUE,
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "text",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans"
  ) #most similar to observed is survivability = 0.11
  venn.diagram(
    x = list(unlist(alliPS), splitBsurvive),
    category.names = c("iPS" , "splitBsim"),
    filename = paste0(plotDirectory,'VD_gDNA/g1MOI0.5/simulation/splitBiPS_survivability',survivability[s],'.pdf'),
    output=TRUE,
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "text",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans"
  ) #most similar to observed is survivability = 0.11
  venn.diagram(
    x = list(splitAsurvive, splitBsurvive),
    category.names = c("splitAsim" , "splitBsim"),
    filename = paste0(plotDirectory,'VD_gDNA/g1MOI0.5/simulation/splitA_splitBoverlap',survivability[s],'.pdf'),
    output=TRUE,
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "text",
    cat.pos = c(-15, 15),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans"
  )
}



#Figure 3C: scRNAseq data
#pull Seurat info
SeuratDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/extractedData/createSeuratobj/"

#no clusters
plusneg_scTransform.integrated <- readRDS(file = paste0(SeuratDirectory,'plusneg_scTransform_50pcs_5000var_filter.rds'))

#clusters res=0.5
plusneg_scTransformclusters_res0.5.integrated = readRDS(file = paste0(SeuratDirectory,'plusneg_scTransform_50pcs_5000var_res0.5_filter.rds'))

###pull out log normalized counts, UMAP coordinates, and Seurat cluster info out of Seurat
plusneglogNormalizedCounts = plusneg_scTransform.integrated[['SCT']]@data #normalized log counts matrix (for counts only, replace "data" by "counts")
plusnegnormalizedCounts = plusneg_scTransform.integrated[['SCT']]@counts #normalized counts matrix, i.e. sequencing depth corrected counts
plusnegcells_count = colnames(plusnegnormalizedCounts) #CellIds with sample name (demultiplexes integrated samples, ex. Wusib1 vs Wusib2) as prefix
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
#mutate(barNum = paste0("B",row_number())) %>% #comment out the above 2 lines and switch with these 2 + saveRDS to save unique barcodes 
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

#make top 100 barcodes table
top100 = plusnegbarcodebreakdown %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  transmute(barcode=barNum, siblingAcells = sibA, siblingBcells = sibB)
top100 = head(top100,100)
print(top100,n=20) #print top 20 barcodes

#color B2, B5, B7, B20
top100_pullouts = c("B2","B5","B7","B20")
top100pullplot = ggplot() +
  geom_point(data=top100,aes(x=siblingBcells,y=siblingAcells))+
  #geom_text(aes(x=siblingBcells+10,label=barcode))+
  geom_point(data=top100[which(top100$barcode %in% top100_pullouts),],aes(x=siblingBcells,y=siblingAcells,color="orange"), size=1)+
  geom_text(data=top100[which(top100$barcode %in% top100_pullouts),],aes(x=siblingBcells,y=siblingAcells+10,label=barcode,color="orange"))+
  xlim(0,310)+ylim(0,310)+
  labs(x="number of cells in sibling B differentiation",
       y="number of cells in sibling A differentiation",
       title="top 100 barcodes associated with the most cells")+
  theme_minimal()+
  theme(legend.position ="none")
#save top 100 barcode scatterplot
pdf(paste(plotDirectory,"top100plot.pdf", sep=""), width=5, height=5)
print(top100pullplot)
dev.off()

#UMAPs: B2, B5, B7, B20
barcodes = c('B2', 'B5', 'B7', 'B20')
for (i in 1: length(barcodes)){
  plusnegbarcode = barcodes[i]
  plusnegbartitle = paste("barcode", plusnegbarcode, ": ")
  barmap = ggplot() +
    geom_point(data = plusnegappendumapCoordinates[which(plusnegappendumapCoordinates$barNum !="NA"),], aes(x = UMAP_1, y = UMAP_2, color = "other"), size = 1) + #all barcoded cells
    geom_point(data = plusnegappendumapCoordinates[which((plusnegappendumapCoordinates$barNum==plusnegbarcode) & (plusnegappendumapCoordinates$SampleNum=="sibA")),], aes(x = UMAP_1, y = UMAP_2, color = "sibA"), size = 3) +
    geom_point(data = plusnegappendumapCoordinates[which((plusnegappendumapCoordinates$barNum==plusnegbarcode) & (plusnegappendumapCoordinates$SampleNum=="sibB")),], aes(x = UMAP_1, y = UMAP_2, color = "sibB"), size = 3) +
    scale_color_manual(breaks = c("sibA","sibB","other"), values = c("sibA"="#FFB6DB","sibB"="#ff068f","other"="gray80"), labels = c(paste("barcode ",plusnegbarcode,"A",sep=""),paste("barcode ",plusnegbarcode,"B",sep=""),"other cells")) +
    create_lpr_theme() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = rel(0.3)),
          legend.text = element_text(size = rel(0.3), angle = 20),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    labs(color =  "Sample origin of barcoded cells") +
    ggtitle(paste0(plusnegbartitle))
  #save pullout UMAPs 
  pdf(paste(plotDirectory,plusnegbarcode,"umap.pdf", sep=""), width=6.5, height=6)
  print(barmap)
  dev.off()
}



