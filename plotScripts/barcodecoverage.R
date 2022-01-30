#Supplemental figure: Barcode coverage
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
plotDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Supplement/coverage/"
plusnegappendumapCoordinates_A = plusnegappendumapCoordinates %>%
  filter(SampleNum == "sibA")

sibAplot <- ggplot() +
  geom_point(data = plusnegappendumapCoordinates_A, aes(x = UMAP_1, y = UMAP_2, color = "other")) +
  geom_point(data = plusnegappendumapCoordinates_A[which((plusnegappendumapCoordinates_A$barNum !="NA")),], aes(x = UMAP_1, y = UMAP_2, color = "sibA"), size = 1) +
  scale_color_manual(breaks = c("sibA","other"), values = c("sibA"="#FFB6DB","other"="gray80"), labels = c("with barcode","without barcode")) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.3)),
        legend.text = element_text(size = rel(0.3), angle = 20),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color =  "Differentiation A cells") +
  ggtitle("where are barcoded cells in differentiation A?") #overlap pretty completely as we would anticipate based on orig.ident distribution

plusnegappendumapCoordinates_B = plusnegappendumapCoordinates %>%
  filter(SampleNum == "sibB")

sibBplot <- ggplot() +
  geom_point(data = plusnegappendumapCoordinates_B, aes(x = UMAP_1, y = UMAP_2, color = "other")) +
  geom_point(data = plusnegappendumapCoordinates_B[which((plusnegappendumapCoordinates_B$barNum !="NA")),], aes(x = UMAP_1, y = UMAP_2, color = "sibB"), size = 1) +
  scale_color_manual(breaks = c("sibB","other"), values = c("sibB"="#ff068f","other"="gray80"), labels = c("with barcode","without barcode")) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.3)),
        legend.text = element_text(size = rel(0.3), angle = 20),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color =  "Differentiation B cells") +
  ggtitle("where are barcoded cells in differentiation B?") #overlap pretty completely as we would anticipate based on orig.ident distribution

#save split A barcode coverage UMAP
pdf(paste(plotDirectory,"barcodesibA_umap.pdf",sep=""), width=6.5, height=6)
print(sibAplot)
dev.off()

#save split B barcode coverage UMAP
pdf(paste(plotDirectory,"barcodesibB_umap.pdf",sep=""), width=6.5, height=6)
print(sibBplot)
dev.off()


#raw barplot
barcodedistributionbysibling <- plusneglogNormalizedCountsSubsetWBarcodes %>%
  group_by(plusnegcells_clusters,SampleNum) %>%
  summarize(barcodedcellsinclusterpersib = n()) %>%
  spread(SampleNum,barcodedcellsinclusterpersib) %>%
  mutate_all(~replace(., is.na(.), 0)) # #make all NAs = 0 (i.e. sibling B cluster 15)

barcodedistributionbysiblingforplot = barcodedistributionbysibling %>% 
  select(plusnegcells_clusters,sibA,sibB) %>%
  gather(SampleNum,barcodedcellsinclusterpersib,sibA:sibB) #%>%


#facet
samplenum.labs <- c(paste0("split A: ",nrow(plusnegappendumapCoordinates_A)," cells"),
                    paste0("split B: ",nrow(plusnegappendumapCoordinates_B)," cells"))
names(samplenum.labs) <- c("sibA", "sibB")

facetplot = ggplot(data=barcodedistributionbysiblingforplot, aes(x=plusnegcells_clusters,y=barcodedcellsinclusterpersib,fill=SampleNum)) +
  geom_bar(stat="identity") +
  facet_wrap(~SampleNum,ncol=1,labeller = labeller(SampleNum = samplenum.labs))+
  ggtitle("All barcoded cells: cluster distribution") +
  ylab("number of cells")+xlab("Seurat clusters, res=0.5") +
  scale_fill_manual(breaks=c("sibA","sibB"), values = c("#FFB6DB","#ff068f")) +
  theme_bw()+
  theme(legend.position = "none")

#save baseline barcoded cell distribution across Seurat clusters
pdf(paste(plotDirectory,"overallclusterdistribution_facet.pdf",sep=""), width=4, height=4)
print(facetplot)
dev.off()



