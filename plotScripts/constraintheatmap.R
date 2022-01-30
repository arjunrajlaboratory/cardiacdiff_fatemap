#Supplemental figure: Heatmap with dendrogram for barcode constraint analysis
#July2021

library(ComplexHeatmap)
library(circlize)

#DEPENDENCY on Figure2_May2021.R: MAKE SURE TO RUN THAT BEFORE CONTINUING
#pull info from Figure 2 directory
dataDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Figure2/"
sibA_20cellbarcodeclusters_new = readRDS(file = paste0(dataDirectory,'sibA_20cellbarcodeclusters_new.rds'))
sibB_20cellbarcodeclusters_new = readRDS(file = paste0(dataDirectory,'sibB_20cellbarcodeclusters_new.rds'))

#plotData directory
plotDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Supplement/constraintheatmap/"

#pull split A observed barcode distributions
p_observed_allA = sibA_20cellbarcodeclusters_new %>% 
  ungroup() %>%
  mutate(prob = prop+0.000001) %>%
  dplyr::select(barNum,plusnegcells_clusters,prob) %>%
  pivot_wider(id_cols = plusnegcells_clusters,names_prefix = "",names_from = barNum,values_from = prob) %>%
  dplyr::select(-plusnegcells_clusters)

p_observed_allA_clust = sibA_20cellbarcodeclusters_new %>% 
  ungroup() %>%
  mutate(prob = prop+0.000001) %>%
  dplyr::select(barNum,plusnegcells_clusters,prob) %>%
  pivot_wider(id_cols = plusnegcells_clusters,names_prefix = "",names_from = barNum,values_from = prob) #keep clusters

p_observed_allA_mat = t(as.matrix(p_observed_allA))
colnames(p_observed_allA_mat) = p_observed_allA_clust$plusnegcells_clusters

#pull split B observed barcode distributions
p_observed_allB = sibB_20cellbarcodeclusters_new %>% 
  ungroup() %>%
  mutate(prob = prop+0.000001) %>%
  dplyr::select(barNum,plusnegcells_clusters,prob) %>%
  pivot_wider(id_cols = plusnegcells_clusters,names_from = barNum,values_from = prob) %>%
  dplyr::select(-plusnegcells_clusters)

p_observed_allB_clust = sibB_20cellbarcodeclusters_new %>% 
  ungroup() %>%
  mutate(prob = prop+0.000001) %>%
  dplyr::select(barNum,plusnegcells_clusters,prob) %>%
  pivot_wider(id_cols = plusnegcells_clusters,names_prefix = "",names_from = barNum,values_from = prob) #keep clusters

p_observed_allB_mat = t(as.matrix(p_observed_allB))
colnames(p_observed_allB_mat) = p_observed_allB_clust$plusnegcells_clusters

#combine datasets into single heatmap
p_observed_all_mat = rbind(p_observed_allA_mat,p_observed_allB_mat)

col_fun = circlize::colorRamp2(c(0, 0.6), c("grey95", "darkblue"))
ha = rowAnnotation(split = c(rep("A",nrow(p_observed_allA_mat)),rep("B",nrow(p_observed_allB_mat))))

combinedheatmap = Heatmap(p_observed_all_mat,name = "normalized cell\nproportion", col = col_fun,
        column_title = "Seurat clusters", column_title_side = "bottom",
        row_title = "all constraint barcodes", row_title_rot = 90,
        cluster_columns = FALSE,
        column_names_rot = 0, 
        clustering_distance_rows = JSD,
        #row_split=4,
        right_annotation = ha) #by default uses euclidean distance, by default agglomerative, complete hierarchical clustering; can set as diana for divisive somehow

#save heatmap showing propensity of each barcode for each Seurat cluster
pdf(paste(plotDirectory,"combinedheatmap.pdf",sep=""), width=6, height=8)
print(combinedheatmap)
dev.off()



