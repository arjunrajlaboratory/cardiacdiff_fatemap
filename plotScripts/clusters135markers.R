#pull markers for Seurat clusters 1/3/5
#11/13/21

library(Seurat)
library(gridExtra)
library(grid)
library(tidyverse)
library(ComplexHeatmap)

#plotData directory
plotDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Supplement/"

#pull Seurat marker info
SeuratDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/extractedData/createSeuratobj/"
#clusters res=0.5
plusneg_scTransformclusters_res0.5.integrated = readRDS(file = paste0(SeuratDirectory,'plusneg_scTransform_50pcs_5000var_res0.5_filter.rds'))

plusneg_scTransformclusters_res0.5.integrated.markers = as_tibble(read.csv(file=paste0(SeuratDirectory,'markers_nonorm.csv')))
cluster1markers = plusneg_scTransformclusters_res0.5.integrated.markers %>% 
  filter(cluster == "1") %>%
  filter(pct.1 >= 0.5) %>%
  arrange(desc(avg_log2FC))
write.csv(file = paste0(plotDirectory,'cluster1markers.csv'), x = cluster1markers)
cluster3markers = plusneg_scTransformclusters_res0.5.integrated.markers %>% 
  filter(cluster == "3") %>%
  filter(pct.1 >= 0.5) %>%
  arrange(desc(avg_log2FC))
write.csv(file = paste0(plotDirectory,'cluster3markers.csv'), x = cluster3markers)
cluster5markers = plusneg_scTransformclusters_res0.5.integrated.markers %>% 
  filter(cluster == "5") %>%
  filter(pct.1 >= 0.5) %>%
  arrange(desc(avg_log2FC))
write.csv(file = paste0(plotDirectory,'cluster5markers.csv'), x = cluster5markers)

#heatmaps
plusneg_scTransformclusters_res0.5.integrated.norm = plusneg_scTransformclusters_res0.5.integrated
DefaultAssay(plusneg_scTransformclusters_res0.5.integrated.norm) <- "RNA" #notice that the zero is not actually at the grey point
NormalizeData(plusneg_scTransformclusters_res0.5.integrated.norm)
cluster_anno <- plusneg_scTransformclusters_res0.5.integrated@meta.data$seurat_clusters
#explicitly map colors to scaled expression values
quantile(mat,c(0.1,0.95))
Seurat::PurpleAndYellow()
## make the black color map to 0. the yellow map to highest and the purle map to the lowest
col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#FF00FF", "black", "#FFFF00"))

#cluster1
mat1 <- plusneg_scTransformclusters_res0.5.integrated.norm[["RNA"]]@data[cluster1markers$gene,] %>% as.matrix()
#after some research decide to use RNA@data slot since that is what is largely
#recommended for differential expression analysis as of 2/12/21: https://github.com/satijalab/seurat/discussions/4032
mat1 <- t(scale(t(mat1)))

#cluster3
mat3 <- plusneg_scTransformclusters_res0.5.integrated.norm[["RNA"]]@data[cluster3markers$gene,] %>% as.matrix()
#after some research decide to use RNA@data slot since that is what is largely
#recommended for differential expression analysis as of 2/12/21: https://github.com/satijalab/seurat/discussions/4032
mat3 <- t(scale(t(mat3)))

#cluster5
mat5 <- plusneg_scTransformclusters_res0.5.integrated.norm[["RNA"]]@data[cluster5markers$gene,] %>% as.matrix()
#after some research decide to use RNA@data slot since that is what is largely
#recommended for differential expression analysis as of 2/12/21: https://github.com/satijalab/seurat/discussions/4032
mat5 <- t(scale(t(mat5)))

#draw heatmap with top 15 markers for each of 3 clusters
top15mat1 = mat1[1:15,]
top15mat3 = mat3[1:15,]
top15mat5 = mat5[1:15,]

top15mat1heatmap_raster4 = Heatmap(top15mat1, name = "Normalized counts",  
                                  column_split = factor(cluster_anno),
                                  cluster_columns = FALSE, #if TRUE will not be in Seurat order
                                  show_column_dend = FALSE,
                                  cluster_column_slices = TRUE,
                                  column_title_gp = gpar(fontsize = 8),
                                  column_gap = unit(0.5, "mm"),
                                  cluster_rows = FALSE,
                                  show_row_dend = TRUE,
                                  col = col_fun,
                                  row_names_gp = gpar(fontsize = 8),
                                  #column_title_rot = 90, #don't rotate 90 degrees
                                  top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(15)))), #match number of clusters
                                  show_column_names = FALSE,
                                  heatmap_legend_param = list(direction = "horizontal"),
                                  use_raster = TRUE,
                                  raster_quality = 4,
                                  row_title = NULL) #increasing raster_quality increases size of file?

top15mat3heatmap_raster4 = Heatmap(top15mat3, name = "Normalized counts",  
                                   column_split = factor(cluster_anno),
                                   cluster_columns = FALSE, #if TRUE will not be in Seurat order
                                   show_column_dend = FALSE,
                                   cluster_column_slices = TRUE,
                                   column_title_gp = gpar(fontsize = 8),
                                   column_gap = unit(0.5, "mm"),
                                   cluster_rows = FALSE,
                                   show_row_dend = TRUE,
                                   col = col_fun,
                                   row_names_gp = gpar(fontsize = 8),
                                   #column_title_rot = 90, #don't rotate 90 degrees
                                   #top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(15)))), #no need bc will be provided by mat1
                                   show_column_names = FALSE,
                                   heatmap_legend_param = list(direction = "horizontal"),
                                   use_raster = TRUE,
                                   raster_quality = 4,
                                   row_title = NULL) #increasing raster_quality increases size of file?

top15mat5heatmap_raster4 = Heatmap(top15mat5, name = "Normalized counts",  
                                   column_split = factor(cluster_anno),
                                   cluster_columns = FALSE, #if TRUE will not be in Seurat order
                                   show_column_dend = FALSE,
                                   cluster_column_slices = TRUE,
                                   column_title_gp = gpar(fontsize = 8),
                                   column_gap = unit(0.5, "mm"),
                                   cluster_rows = FALSE,
                                   show_row_dend = TRUE,
                                   col = col_fun,
                                   row_names_gp = gpar(fontsize = 8),
                                   #column_title_rot = 90, #don't rotate 90 degrees
                                   #top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(15)))), #no need bc will be provided by mat1
                                   show_column_names = FALSE,
                                   heatmap_legend_param = list(direction = "horizontal"),
                                   use_raster = TRUE,
                                   raster_quality = 4,
                                   row_title = NULL) #increasing raster_quality increases size of file?

htlist = top15mat1heatmap_raster4 %v% top15mat3heatmap_raster4 %v% top15mat5heatmap_raster4

#save htlist
pdf(paste0(plotDirectory,"htlist_raster4.pdf"), width=6, height=6)
draw(htlist, heatmap_legend_side = "bottom",row_gap = unit(0.1, "mm"),row_title = NULL)
dev.off()