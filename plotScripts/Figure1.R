#Figure 1: Seurat clusters, heat map, UMAP, latent dimension
#May2021

library(Seurat)
library(gridExtra)
library(grid)
library(tidyverse)
library(ComplexHeatmap)
library(presto) # https://github.com/immunogenomics/presto
library(tictoc)
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

#plotData directory
plotDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Figure1/"

#Figure 1B: Seurat clusters, res=0.5
#calculate centroid to label clusters
addcentroids_clusters = plusnegumapCoordinates %>%
  group_by(plusnegcells_clusters) %>%
  summarize(UMAP_1 = mean(UMAP_1),UMAP_2=mean(UMAP_2)) 

Seuratclust = ggplot() +
  geom_point(plusnegumapCoordinates,mapping = aes(x = UMAP_1, y = UMAP_2,color = plusnegcells_clusters),size = 1) +
  geom_text(mapping = aes_string(x = addcentroids_clusters$UMAP_1, 
                                 y = addcentroids_clusters$UMAP_2),
            color = "black", size = 8,
            label=addcentroids_clusters$plusnegcells_clusters)+
  create_lpr_theme() + #function written by Lee (see bottom)
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color =  "Seurat-defined\n clusters") +
  ggtitle("Seurat-defined clusters, res=0.5")
#save differentiated cell expression UMAP colored by Seurat cluster
pdf(paste0(plotDirectory,"Figure1B.pdf"), width=6, height=5)
print(Seuratclust)
dev.off()

#Figure 1C: 
#Heatmap of markers: run first 2 commands to create marker list then load afterwards so as not to have to wait
plusneg_scTransformclusters_res0.5.integrated.markers <- FindAllMarkers(plusneg_scTransformclusters_res0.5.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
write.csv(file = paste0(SeuratDirectory,'markers_nonorm.csv'), x = plusneg_scTransformclusters_res0.5.integrated.markers)
plusneg_scTransformclusters_res0.5.integrated.markers = as_tibble(read.csv(file=paste0(SeuratDirectory,'markers_nonorm.csv')))
selectedx = c('DCN','LUM','DDR2','MXRA8','SEMA6A','HEY1','TBX1','CALM1','TNNT2','ACTC1',
              'MYL7','CRB21','MT-ND21','RRM11','SCLT1','ACAT2','PCLAF','PLCB11','ISL1',
              'TBX18','WT1','KPNA2','MKI671','HES42','ANXA1','EPCAM','CDH1',
              'NR2F1','LAMP5','CSKMT1','PLCG2')
#pick 2-3 per cluster:
#0:'DCN','LUM'
#1:'DDR2','MXRA8'
#2:'SEMA6A','HEY1'
#3:'TBX1','CALM1'
#4:'TNNT2','ACTC1','MYL7'
#5:'CRB21','MT-ND21'
#6:'RRM11','SCLT1'
#7:'ACAT2','PCLAF'
#8:'PLCB11','ISL1'
#9:'TBX18','WT1'
#10:'KPNA2','MKI671'
#11:'HES42','ANXA1'
#12:'EPCAM','CDH1'
#13:'NR2F1','LAMP5'
#14:'CSKMT1','PLCG2'
select2 <- plusneg_scTransformclusters_res0.5.integrated.markers %>%
  filter(X %in% selectedx) %>%
  filter(avg_log2FC>1)

#initially used Seurat command but then switched to ComplexHeatmap to allow for clustering info
#based on Ming Tang's post: https://divingintogeneticsandgenomics.rbind.io/post/enhancement-of-scrnaseq-heatmap-using-complexheatmap/
plusneg_scTransformclusters_res0.5.integrated.norm = plusneg_scTransformclusters_res0.5.integrated
DefaultAssay(plusneg_scTransformclusters_res0.5.integrated.norm) <- "RNA" #notice that the zero is not actually at the grey point
NormalizeData(plusneg_scTransformclusters_res0.5.integrated.norm)
mat <- plusneg_scTransformclusters_res0.5.integrated.norm[["RNA"]]@data[select2$gene,] %>% as.matrix()
#after some research decide to use RNA@data slot since that is what is largely
  #recommended for differential expression analysis as of 2/12/21: https://github.com/satijalab/seurat/discussions/4032
mat <- t(scale(t(mat)))

cluster_anno <- plusneg_scTransformclusters_res0.5.integrated@meta.data$seurat_clusters

#explicitly map colors to scaled expression values
quantile(mat,c(0.1,0.95))
Seurat::PurpleAndYellow()
## make the black color map to 0. the yellow map to highest and the purle map to the lowest
col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#FF00FF", "black", "#FFFF00"))
markertype = c(rep(0,2),rep(1,2),rep(2,2),rep(3,2),rep(4,3),
               rep(5,2),rep(6,2),rep(7,2),rep(8,2),rep(9,2),
               rep(10,2),rep(11,2),rep(12,2),rep(13,2),rep(14,2))

selectmarkerheatmap_raster4 = Heatmap(mat, name = "Normalized counts",  
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

#save heatmap demonstrating expression of select markers in differentiated cells
pdf(paste0(plotDirectory,"Figure1Cheatmap_raster4.pdf"), width=6, height=6)
draw(selectmarkerheatmap_raster4, heatmap_legend_side = "bottom",row_gap = unit(0.1, "mm"),
     row_split = factor(markertype),row_title = NULL)
dev.off()

#UMAPs: walk through list of genes: TNNT2, ISL1, EPCAM, LUM, WT1
genes = c('TNNT2','ISL1','EPCAM','LUM','WT1')
for (i in 1:length(genes)){
  gene = genes[i]
  geneumap = ggplot(plusnegumapCoordinates, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes_string(color = pull(plusneglogNormalizedCounts,genes[i])), size = 1) +
    scale_color_gradient(low = "gray95", high = "darkblue") +
    create_lpr_theme() + #function written by Lee (see bottom)
    theme(legend.position = "bottom",
          legend.title = element_text(size = rel(0.6)),
          legend.text = element_text(size = rel(0.6), angle = 30),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    labs(color =  "Scaled\n Log UMI counts") +
    ggtitle(paste(genes[i],"expression"))
  #save UMAPs displaying marker expression
  pdf(paste(plotDirectory,genes[i],"_umap.pdf",sep=""), width=6.5, height=6)
  print(geneumap)
  dev.off()
}

#Figure 1D
#See LatentDimension...R document for more info
#iPSC
# Load the relevant datasets (filteredbarcode matrix files)
iPSCBalone <- readRDS(file = paste0(SeuratDirectory,'iPSCBaloneSCT_5000var_filter.rds'))
iPSCB_PCA <- readRDS(file = paste0(SeuratDirectory,'iPSCB_PCA_5000var_filter.rds'))

#follow zeehio's response on https://github.com/satijalab/seurat/issues/982/
mat_iPS <- GetAssayData(iPSCBalone,assay='SCT',slot="scale.data") #rows are genes, columns are cells
pca_iPS <- iPSCB_PCA[['pca']]
total_variance_iPS <- sum(matrixStats::rowVars(mat_iPS)) #find variance estimate per row, i.e. per gene
#4581.422

eigValues_iPS=(pca_iPS@stdev)^2
varExplained_iPS = as_tibble(eigValues_iPS/total_variance_iPS)%>% #sum is 0.28 (<<1 bc this is only first 50 PCs maybe??)
  transmute(varExplained = value) %>%
  mutate(PC = row_number(),cumulative = cumsum(varExplained))

iPSvar = ggplot(varExplained_iPS,aes(x=PC))+
  geom_bar(aes(y=varExplained),stat='identity')+
  ylim(0,0.1)+
  xlab("PC")+
  ylab("fraction of variance explained")+
  ggtitle("iPS sample variance explained")+
  theme_bw()

#randomize by mixing up scale_data from iPS
mat_iPSrandom <- mat_iPS #ty to Lauren Beck for help
set.seed(2059)
for (index_row in 1:nrow(mat_iPSrandom)){
  mat_iPSrandom[index_row,] <- sample(mat_iPSrandom[index_row,], ncol(mat_iPS))
}
iPSCBrandom <- SetAssayData(
  object = iPSCBalone,
  slot = "scale.data",
  new.data = mat_iPSrandom,
  assay = "SCT"
) 
iPSCBrandom <- RunPCA(iPSCBrandom, assay='SCT',features = VariableFeatures(object = iPSCBrandom), verbose = FALSE)

pca_iPSrandom <- iPSCBrandom[['pca']]
total_variance_iPSrandom <- sum(matrixStats::rowVars(mat_iPSrandom)) #find variance estimate per row, i.e. per gene

eigValues_iPSrandom=(pca_iPSrandom@stdev)^2
varExplained_iPSrandom = as_tibble(eigValues_iPSrandom/total_variance_iPSrandom)%>% #sum is 0.729 (<1 bc this is only first 50 PCs maybe)
  transmute(varExplained = value) %>%
  mutate(PC = row_number(),cumulative = cumsum(varExplained))

iPSvarrandom = ggplot(varExplained_iPSrandom,aes(x=PC))+
  geom_bar(aes(y=varExplained),stat='identity')+
  ylim(0,0.3)+
  xlab("PC")+
  ylab("fraction of variance explained")+
  ggtitle("iPS randomized variance explained")+
  theme_bw()
#add to original iPS plot
iPSvar_line = iPSvar +
  geom_line(varExplained_iPSrandom,mapping=aes(x=PC,y=varExplained),stat='identity',color='#F8766D')+
  xlim(0,10)

#diff: use integrated data that we have been using
mat_diff <- GetAssayData(plusneg_scTransform.integrated,assay='SCT',slot="scale.data")
pca_diff <- plusneg_scTransform.integrated[['pca']]
total_variance_diff <- sum(matrixStats::rowVars(mat_diff)) #4926.843

eigValues_diff=(pca_diff@stdev)^2
varExplained_diff = as_tibble(eigValues_diff/total_variance_diff)%>% #sum is 0.261 (<<1 bc this is only first 50 PCs maybe)
  transmute(varExplained = value) %>%
  mutate(PC = row_number(),cumulative = cumsum(varExplained))

diffvar = ggplot(varExplained_diff,aes(x=PC))+
  geom_bar(aes(y=varExplained),stat='identity')+
  ylim(0,0.1)+
  xlab("PC")+
  ylab("fraction of variance explained")+
  ggtitle("differentiated sample variance explained")+
  theme_bw()

#randomize as we did for iPS
mat_diffrandom <- mat_diff 
set.seed(2059)
for (index_row in 1:nrow(mat_diffrandom)){
  mat_diffrandom[index_row,] <- sample(mat_diffrandom[index_row,], ncol(mat_diff))
}
diffrandom <- SetAssayData(
  object = plusneg_scTransform.integrated,
  slot = "scale.data",
  new.data = mat_diffrandom,
  assay = "SCT"
) 
diffrandom <- RunPCA(diffrandom, assay='SCT',features = VariableFeatures(object = diffrandom), verbose = FALSE)

mat_diffrandom <- as.matrix(GetAssayData(diffrandom,assay='SCT',slot="scale.data")) #rows are genes, columns are cells
#same as mat_diff_SD_random
pca_diffrandom <- diffrandom[['pca']]
total_variance_diffrandom <- sum(matrixStats::rowVars(mat_diffrandom)) #find variance estimate per row, i.e. per gene

eigValues_diffrandom=(pca_diffrandom@stdev)^2
varExplained_diffrandom = as_tibble(eigValues_diffrandom/total_variance_diffrandom)%>% #sum is 0.729 (<1 bc this is only first 50 PCs maybe)
  transmute(varExplained = value) %>%
  mutate(PC = row_number(),cumulative = cumsum(varExplained))

diffvarrandom = ggplot(varExplained_diffrandom,aes(x=PC))+
  geom_bar(aes(y=varExplained),stat='identity')+
  ylim(0,0.3)+
  xlab("PC")+
  ylab("fraction of variance explained")+
  ggtitle("diff randomized variance explained")+
  theme_bw()
#add to original diff plot
diffvar_line = diffvar +
  geom_line(varExplained_diffrandom,mapping=aes(x=PC,y=varExplained),stat='identity',color='#F8766D')+
  xlim(0,10)

#save plot comparing fraction of variance explained by top 10 hiPS cell PCs
#vs fraction of variance explained by top 10 differentiated cell PCs
pdf(paste0(plotDirectory,"Figure1D.pdf"), width=6, height=6)
print(grid.arrange(iPSvar_line,diffvar_line,nrow=2))
dev.off()


