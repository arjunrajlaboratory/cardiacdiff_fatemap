#Run all plots
#July 2021

#make sure to run extraction scripts before running plot scripts
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/extractionScripts/createSeuratobject.R")
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/extractionScripts/createbarcodeobject.R")
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/extractionScripts/gDNAbarcodeobjects.R")

#Figure 1
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/Figure1.R")

#Figure 2: must be run before Supplementary Figure 2
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/Figure2.R")

#Figure 3: must be run before Supplementary Figure 3
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/Figure3.R")

#Figure 4:
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/Figure4.R")

#Figure 5:
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/Figure5.R")


#Supplement
#Supp Figure 1:
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/additionalMarkers.R")

#Supp Figure 2:
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/GFPdistribution.R")

#Supp Figure 3: 
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/barcodecoverage.R")
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/proportionbarcodedpercluster.R")
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/clusters135.R")

#Supp Figure 4:
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/barcodeprobabilitydistributiongeneration.R")
#below script requires RDS files created in Figure 2 script
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/constraintheatmap.R")
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/cluster3significance.R")

#Supp Figure 5:
#below script requires RDS files created in Figure 3 script
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/survivaloverlapsignificance.R")
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/survivalVD.R")
print(top100,n=20)

#Supp Figure 6:
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/Figure4Bchangecutoffs.R")
source("/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotScripts/JSDmelanomadata.R")
