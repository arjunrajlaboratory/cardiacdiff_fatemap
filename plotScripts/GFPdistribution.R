#GFP positive vs negative distribution
#10/31/21

#"GFP positive and negative populations are indistinguishable"
library(tidyverse)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(grid)
library(DescTools)

#grab plusnegappendumapCoordinates from all Seurat scripts
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

#plotData directory
plotDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Supplement/"

#umap
plusnegappendumapCoordinates_GFP = plusnegappendumapCoordinates %>%
  filter(SampleNum == "sibA"|SampleNum == "sibB") #10961 GFP+ of 17599 total profiled cells

GFPplot <- ggplot() +
  geom_point(data = plusnegappendumapCoordinates, aes(x = UMAP_1, y = UMAP_2, color = "other"), size = 1) +
  geom_point(data = plusnegappendumapCoordinates_GFP, aes(x = UMAP_1, y = UMAP_2, color = "GFP"), size = 1) +
  scale_color_manual(breaks = c("GFP","other"), values = c("GFP"="darkgreen","other"="gray80"), labels = c("GFP positive","GFP negative")) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.3)),
        legend.text = element_text(size = rel(0.3), angle = 20),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color =  "Sequenced cells") +
  ggtitle("GFP positive vs GFP negative cell distribution") #overlap pretty closely, i.e. "indistinguishable"

pdf(paste(plotDirectory,"GFPindistinguishable_umap.pdf",sep=""), width=6.5, height=6)
print(GFPplot)
dev.off()

#bar graph
GFPdistribution <- plusneglogNormalizedCounts %>%
  mutate(GFP = case_when(
    SampleNum=="negative" ~ "GFPneg",
    SampleNum=="sibA" | SampleNum=="sibB" ~ "GFPpos"
  )) %>%
  group_by(plusnegcells_clusters,GFP) %>%
  summarize(cellsinclusterperGFP = n()) %>%
  spread(GFP,cellsinclusterperGFP) %>%
  mutate_all(~replace(., is.na(.), 0)) # #make all NAs = 0 (i.e. sibling B cluster 15)

GFPdistributionforplot = GFPdistribution %>% 
  select(plusnegcells_clusters,GFPneg,GFPpos) %>%
  gather(GFP,cellsinclusterperGFP,GFPneg:GFPpos) #%>%


#facet
GFPnum.labs <- c(paste0("GFP negative: ",nrow(plusnegappendumapCoordinates)-nrow(plusnegappendumapCoordinates_GFP)," cells"),
                    paste0("GFP positive: ",nrow(plusnegappendumapCoordinates_GFP)," cells"))
names(GFPnum.labs) <- c("GFPneg", "GFPpos")

GFPfacetplot = ggplot(data=GFPdistributionforplot, aes(x=plusnegcells_clusters,y=cellsinclusterperGFP,fill=GFP)) +
  geom_bar(stat="identity") +
  facet_wrap(~GFP,ncol=1,labeller = labeller(GFP = GFPnum.labs))+
  ggtitle("All cells: GFP distribution") +
  ylab("number of cells")+xlab("Seurat clusters, res=0.5") +
  scale_fill_manual(breaks=c("GFPneg","GFPpos"), values = c("grey80","darkgreen")) +
  theme_bw()+
  theme(legend.position = "none")

#save baseline barcoded cell distribution across Seurat clusters
pdf(paste(plotDirectory,"GFPclusterdistribution_facet.pdf",sep=""), width=4, height=4)
print(GFPfacetplot)
dev.off()



