#Seurat extraction script for paper
#cycle through resolutions
#last edited July 2021

library(Seurat)
library(grid)
library(gridExtra)

dataDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/extractedData/createSeuratobj/"

# create Seurat object for differentiated cells
# Load the relevant datasets (filteredbarcode matrix files)
siblingA.data <- Read10X(data.dir = "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/rawData/10XData/cellranger_outs/20210224_10kbarsiblingA/filtered_feature_bc_matrix")
siblingB.data <- Read10X(data.dir = "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/rawData/10XData/cellranger_outs/20210224_10kbarsiblingB/filtered_feature_bc_matrix")
negative.data <- Read10X(data.dir = "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/rawData/10XData/cellranger_outs/20210224_10kbarsiblingBneg/filtered_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data).
siblingA <- CreateSeuratObject(counts = siblingA.data, project = "siblingA", min.cells = 3, min.features = 200)
siblingB <- CreateSeuratObject(counts = siblingB.data, project = "siblingB", min.cells = 3, min.features = 200)
negative <- CreateSeuratObject(counts = negative.data, project = "negative", min.cells = 3, min.features = 200)

#help to choose nFeature cutoffs and see % mito and ribo
siblingA[["percent.mt"]] <- PercentageFeatureSet(siblingA, pattern = "^MT-")
siblingA[["percent.ribo"]] <- PercentageFeatureSet(siblingA, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
VlnPlot(siblingA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
siblingB[["percent.mt"]] <- PercentageFeatureSet(siblingB, pattern = "^MT-")
siblingB[["percent.ribo"]] <- PercentageFeatureSet(siblingB, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
VlnPlot(siblingB, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
negative[["percent.mt"]] <- PercentageFeatureSet(negative, pattern = "^MT-")
negative[["percent.ribo"]] <- PercentageFeatureSet(negative, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
VlnPlot(negative, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)

siblingA <- subset(siblingA, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
siblingB <- subset(siblingB, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
negative <- subset(negative, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)

plusneg_scTransform <- merge(x = negative, y = c(siblingA, siblingB), add.cell.ids = c("negative","sibA","sibB"), project = "negsibAsibB") #17600 cells with 24728 features

plusneg_scTransform.list <- SplitObject(plusneg_scTransform, split.by = "orig.ident")
for (i in 1:length(plusneg_scTransform.list)) {
  plusneg_scTransform.list[[i]] <- SCTransform(plusneg_scTransform.list[[i]], verbose = FALSE)
}

options(future.globals.maxSize = 4000 * 1024^2) #increased future.globals.maxSize to 4gb to allow PrepSCT, must be rerun whenever R restarts

varfeatures = c(1000,2000,5000,7000) #7000 is on the high end since there aren't that many cells with 7500+ genes reported
for (j in 1:length(varfeatures)){
  plusneg_scTransform.features <- SelectIntegrationFeatures(object.list = plusneg_scTransform.list, nfeatures = varfeatures[j])
  plusneg_scTransform.list <- PrepSCTIntegration(object.list = plusneg_scTransform.list,
                                                 anchor.features = plusneg_scTransform.features,
                                                 verbose = FALSE)
  plusneg_scTransform.anchors <- FindIntegrationAnchors(object.list = plusneg_scTransform.list,
                                                        normalization.method = "SCT",
                                                        anchor.features = plusneg_scTransform.features,
                                                        verbose = FALSE)
  plusneg_scTransform.integrated <- IntegrateData(anchorset = plusneg_scTransform.anchors,
                                                  normalization.method = "SCT",
                                                  verbose = FALSE)
  plusneg_scTransform.integrated <- RunPCA(plusneg_scTransform.integrated, verbose = FALSE)
  plusneg_scTransform.integrated <- RunUMAP(plusneg_scTransform.integrated, dims = 1:50) #enter this below into logNormalizedCounts
  print(table(Idents(plusneg_scTransform.integrated))) #shows breakdown by sample number of how many cells
  #sibA has 6112, sibB has 4849, neg when included has 6639 cells
  varumap = DimPlot(plusneg_scTransform.integrated)
  print(varumap)
  saveRDS(plusneg_scTransform.integrated, file = paste0(dataDirectory,'plusneg_scTransform_50pcs_',varfeatures[j],'var_filter.rds'))
  }
#look pretty similar to each other; almost no difference between 2000-7000 variable features

#vary resolutions on with same number of variable features: use 5000 variables since it looks basically identical to 7000 var
plusneg_scTransform.integrated = readRDS(paste0(dataDirectory,'plusneg_scTransform_50pcs_5000var_filter.rds'))

res = c(0.4,0.5,0.6,0.8,1.0,1.2)
for (k in 1:length(res)){
  plusneg_scTransformclusters.integrated <- FindNeighbors(plusneg_scTransform.integrated, reduction = "pca", dims = 1:50)
  plusneg_scTransformclusters.integrated <- FindClusters(plusneg_scTransformclusters.integrated, resolution = res[k]) 
  resplot <- DimPlot(plusneg_scTransformclusters.integrated, label = TRUE, label.size=10)
  print(resplot)
  saveRDS(plusneg_scTransformclusters.integrated, file = paste0(dataDirectory,'plusneg_scTransform_50pcs_5000var_res',res[k],'_filter.rds'))
  }

TNNT2 = FeaturePlot(plusneg_scTransformclusters.integrated, features="TNNT2")
ACTC1 = FeaturePlot(plusneg_scTransformclusters.integrated, features="ACTC1")
MKI67 = FeaturePlot(plusneg_scTransformclusters.integrated, features="MKI67")

#use 5000 variable features and resolution = 0.5 for analyses based on TNNT2/ACTC1 cluster and MKI67 cluster

# create Seurat object for iPS cells
# Load the relevant datasets (filteredbarcode matrix files)
iPSCB.data <- Read10X(data.dir = "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/rawData/10XData/cellranger_outs/20200810_iPSCsample8/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
iPSCB <- CreateSeuratObject(counts = iPSCB.data, project = "iPSCB", min.cells = 3, min.features = 200)
iPSCBsub <- subset(iPSCB, subset = nFeature_RNA > 200 & nFeature_RNA < 10000) #consider subsetting out by mt
#supposedly SCTransform replaces NormalizeData, FindVariableFeatures, and ScaleData
iPSCBalone <- SCTransform(object = iPSCBsub, variable.features.n = 5000, verbose = FALSE)
iPSCB_PCA <- RunPCA(iPSCBalone, assay='SCT',features = VariableFeatures(object = iPSCBalone), verbose = FALSE)
saveRDS(iPSCBalone, file = paste0(dataDirectory,'iPSCBaloneSCT_5000var_filter.rds'))
saveRDS(iPSCB_PCA, file = paste0(dataDirectory,'iPSCB_PCA_5000var_filter.rds'))


